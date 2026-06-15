#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import sys
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


EARLY_STAGES = ("I", "II")
LATE_STAGES = ("III", "IV", "III or IV")
DEFAULT_EPITHELIAL_TYPES = (
    "Alveolar cell type 1",
    "Alveolar cell type 2",
    "ROS1+ healthy epithelial",
    "transitional club/AT2",
    "Club",
    "Ciliated",
)


def read_structured(path: str | pathlib.Path) -> dict[str, Any]:
    path = pathlib.Path(path)
    with path.open() as handle:
        if path.suffix.lower() in {".yaml", ".yml"}:
            return yaml.safe_load(handle)
        return json.load(handle)


def write_json(path: pathlib.Path, value: Any) -> None:
    path.write_text(json.dumps(as_jsonable(value), indent=2, sort_keys=True))


def compatible_dataframe(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result.index = pd.Index(
        np.asarray(result.index.astype(str), dtype=object),
        dtype=object,
    )
    for column in result.columns:
        values = result[column]
        if (
            isinstance(values.dtype, pd.StringDtype)
            or isinstance(values.dtype, pd.CategoricalDtype)
            or pd.api.types.is_object_dtype(values.dtype)
        ):
            result[column] = pd.Series(
                [None if pd.isna(value) else str(value) for value in values],
                index=result.index,
                dtype=object,
            )
    return result


def write_compatible_h5ad(adata: ad.AnnData, path: pathlib.Path) -> None:
    result = adata.copy()
    result.obs = compatible_dataframe(result.obs)
    result.var = compatible_dataframe(result.var)
    result.write_h5ad(path)


def as_jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, pathlib.Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): as_jsonable(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [as_jsonable(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return str(value)


def package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def canonicalize_stage(value: Any, replacements: Mapping[str, str] | None = None) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "none", "unknown", "not available"}:
        return None
    if replacements and text in replacements:
        text = str(replacements[text]).strip()
    normalized = re.sub(r"[\s_-]+", " ", text.upper()).strip()
    if normalized in {"III OR IV", "III/IV", "III IV"}:
        return "III or IV"
    match = re.match(r"^(IV|III|II|I)(?:[ABC]|\d|$)", normalized)
    return match.group(1) if match else None


def normalize_stage_series(
    values: pd.Series,
    replacements: Mapping[str, str] | None = None,
) -> pd.Series:
    return values.map(lambda value: canonicalize_stage(value, replacements)).astype("object")


def select_subcluster_cells(
    obs: pd.DataFrame,
    *,
    cell_type_key: str = "cell_type",
    tumor_pattern: str = "Tumor",
    epithelial_types: Sequence[str] = DEFAULT_EPITHELIAL_TYPES,
) -> pd.Series:
    if cell_type_key not in obs:
        raise KeyError(f"Missing required cell type column {cell_type_key!r}")
    labels = obs[cell_type_key].astype("object").fillna("").astype(str)
    return labels.str.contains(tumor_pattern, regex=False) | labels.isin(epithelial_types)


def split_stage_masks(obs: pd.DataFrame, stage_key: str = "stage") -> dict[str, pd.Series]:
    if stage_key not in obs:
        raise KeyError(f"Missing required stage column {stage_key!r}")
    stages = obs[stage_key].astype("object")
    return {
        "early": stages.isin(EARLY_STAGES),
        "late": stages.isin(LATE_STAGES),
        "unknown": ~stages.isin(EARLY_STAGES + LATE_STAGES),
    }


def effective_faiss_params(
    n_obs: int,
    n_neighbors: int,
    nlist: int | None = None,
    nprobe: int | None = None,
) -> tuple[int, int, int]:
    if n_obs < 3:
        raise ValueError("Subclustering requires at least three cells")
    effective_neighbors = min(max(2, int(n_neighbors)), n_obs - 1)
    effective_nlist = int(nlist) if nlist else min(4096, max(1, round(math.sqrt(n_obs))))
    effective_nlist = min(max(1, effective_nlist), n_obs)
    effective_nprobe = int(nprobe) if nprobe else min(64, effective_nlist)
    effective_nprobe = min(max(1, effective_nprobe), effective_nlist)
    return effective_neighbors, effective_nlist, effective_nprobe


def symmetric_affinity_graph(
    indices: np.ndarray,
    distances: np.ndarray,
    n_obs: int,
) -> sparse.csr_matrix:
    if indices.shape != distances.shape:
        raise ValueError("Neighbor indices and distances must have matching shapes")
    valid = (indices >= 0) & (indices < n_obs)
    positive = distances[valid & (distances > 0)]
    scale = float(np.median(positive)) if positive.size else 1.0
    scale = max(scale, np.finfo(np.float32).eps)
    rows = np.repeat(np.arange(n_obs), indices.shape[1])[valid.ravel()]
    cols = indices.ravel()[valid.ravel()]
    weights = np.exp(-distances.ravel()[valid.ravel()] / scale).astype(np.float32)
    graph = sparse.coo_matrix((weights, (rows, cols)), shape=(n_obs, n_obs)).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    return graph.maximum(graph.T).tocsr()


def validate_mapping(labels: Sequence[Any], mapping: Mapping[str, str]) -> dict[str, str]:
    normalized = {str(key): str(value) for key, value in mapping.items()}
    observed = {str(value) for value in labels}
    missing = sorted(observed - set(normalized))
    extra = sorted(set(normalized) - observed)
    if missing or extra:
        raise ValueError(f"Mapping mismatch; missing clusters={missing}, extra clusters={extra}")
    return normalized


def apply_mapping_exports(
    merged: ad.AnnData,
    stage_name: str,
    stage_indices: Sequence[str],
    cluster_labels: Sequence[Any],
    mapping: Mapping[str, str],
    output_dir: pathlib.Path,
) -> list[pathlib.Path]:
    mapping = validate_mapping(cluster_labels, mapping)
    adjusted = merged.obs["cell_type"].astype(str).copy()
    mapped = pd.Series(
        [mapping[str(label)] for label in cluster_labels],
        index=pd.Index(stage_indices),
        dtype="object",
    )
    adjusted.loc[mapped.index] = mapped
    written: list[pathlib.Path] = []
    for dataset, dataset_obs in merged.obs.groupby("dataset", sort=True, observed=True):
        frame = pd.DataFrame(
            {
                "index": dataset_obs["source_index"].astype(str),
                "cell_type_adjusted": adjusted.loc[dataset_obs.index].astype(str).to_numpy(),
                "cell_type": dataset_obs["cell_type"].astype(str).to_numpy(),
                "batch": dataset_obs["batch"].astype(str).to_numpy(),
            }
        )
        path = output_dir / f"{dataset}_predicted_leiden_{stage_name}.csv"
        frame.to_csv(path, index=False)
        written.append(path)
    return written


def _dense_float32(value: Any) -> np.ndarray:
    if sparse.issparse(value):
        value = value.toarray()
    return np.asarray(value, dtype=np.float32)


def _read_obs_csv(entry: Mapping[str, Any], expected_rows: int, latent_index: pd.Index) -> pd.DataFrame | None:
    obs_csv = entry.get("obs_csv")
    if not obs_csv:
        return None
    index_col = entry.get("obs_csv_index_col", 0)
    frame = pd.read_csv(obs_csv, index_col=index_col)
    join_mode = entry.get("obs_join", "index")
    if join_mode == "position":
        if len(frame) != expected_rows:
            raise ValueError(f"{obs_csv} has {len(frame)} rows; expected {expected_rows}")
        frame = frame.copy()
        frame.index = latent_index
        return frame
    if join_mode != "index":
        raise ValueError(f"Unsupported obs_join {join_mode!r}; use index or position")
    missing = latent_index.difference(frame.index)
    if len(missing):
        raise ValueError(f"{obs_csv} is missing {len(missing)} latent cell IDs")
    return frame.loc[latent_index].copy()


def _resolve_column(
    frames: Sequence[pd.DataFrame],
    key: str | None,
    fallbacks: Sequence[str],
    default: Any = None,
) -> pd.Series:
    for frame in frames:
        for candidate in [key, *fallbacks]:
            if candidate and candidate in frame:
                return frame[candidate].copy()
    if default is not None:
        return pd.Series([default] * len(frames[0]), index=frames[0].index, dtype="object")
    raise KeyError(f"None of the required columns exist: {[key, *fallbacks]}")


def load_manifest_entry(entry: Mapping[str, Any]) -> ad.AnnData:
    name = str(entry["name"])
    latent = ad.read_h5ad(entry["latent_h5ad"], backed="r")
    representation = entry.get("representation", "X")
    if representation == "X":
        matrix = _dense_float32(latent.X)
    elif representation in latent.obsm:
        matrix = _dense_float32(latent.obsm[representation])
    else:
        raise KeyError(f"{entry['latent_h5ad']} has no representation {representation!r}")

    latent_obs = latent.obs.copy()
    external_obs = _read_obs_csv(entry, latent.n_obs, latent.obs_names)
    frames = [frame for frame in [external_obs, latent_obs] if frame is not None]
    cell_type = _resolve_column(
        frames,
        entry.get("cell_type_key"),
        ("predicted_cell_type", "cell_type", "cell_type_tumor", "predictions"),
    )
    batch = _resolve_column(
        frames,
        entry.get("batch_key"),
        ("sample", "batch", "dataset"),
        default=name,
    )
    if entry.get("stage_override") is not None:
        stage_raw = pd.Series([entry["stage_override"]] * latent.n_obs, index=latent.obs_names)
    else:
        stage_raw = _resolve_column(
            frames,
            entry.get("stage_key"),
            ("stage", "uicc_stage", "Pathological stage", "Clinical Stage"),
            default=None,
        )
    stage = normalize_stage_series(stage_raw, entry.get("stage_replacements"))

    source_index = pd.Index(latent.obs_names.astype(str))
    if external_obs is not None and entry.get("source_index_key") in external_obs:
        source_index = pd.Index(external_obs[entry["source_index_key"]].astype(str))
    internal_index = pd.Index([f"{name}::{value}" for value in source_index])
    obs = pd.DataFrame(
        {
            "source_index": source_index.to_numpy(),
            "cell_type": cell_type.astype(str).to_numpy(),
            "dataset": name,
            "batch": batch.astype(str).to_numpy(),
            "stage": stage.to_numpy(),
        },
        index=internal_index,
    )
    latent.file.close()
    return ad.AnnData(X=matrix, obs=obs)


def validate_premerged(adata: ad.AnnData) -> ad.AnnData:
    required = {"cell_type", "dataset", "batch", "stage"}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise ValueError(f"Premerged input is missing obs columns: {missing}")
    result = ad.AnnData(X=_dense_float32(adata.X), obs=adata.obs.copy())
    if "source_index" not in result.obs:
        result.obs["source_index"] = result.obs_names.astype(str)
    result.obs["stage"] = normalize_stage_series(result.obs["stage"])
    result.obs_names = pd.Index(
        [f"{dataset}::{source}" for dataset, source in zip(result.obs["dataset"], result.obs["source_index"])]
    )
    if not result.obs_names.is_unique:
        result.obs_names_make_unique()
    return result


def load_cohort(cfg: Mapping[str, Any]) -> ad.AnnData:
    input_h5ad = cfg.get("subcluster_input_h5ad")
    manifest_path = cfg.get("subcluster_cohort_manifest")
    if bool(input_h5ad) == bool(manifest_path):
        raise ValueError("Set exactly one of subcluster_input_h5ad or subcluster_cohort_manifest")
    if input_h5ad:
        return validate_premerged(ad.read_h5ad(input_h5ad))

    manifest_path = pathlib.Path(manifest_path).resolve()
    manifest = read_structured(manifest_path)
    entries = []
    if manifest.get("atlas"):
        entries.append(manifest["atlas"])
    entries.extend(manifest.get("queries") or [])
    if not entries:
        raise ValueError("Cohort manifest contains no atlas or query entries")
    resolved_entries = []
    for entry in entries:
        entry = dict(entry)
        for key in ("latent_h5ad", "obs_csv"):
            if entry.get(key):
                path = pathlib.Path(entry[key])
                entry[key] = str(path if path.is_absolute() else manifest_path.parent / path)
        resolved_entries.append(entry)
    cohorts = [load_manifest_entry(entry) for entry in resolved_entries]
    dimensions = {cohort.n_vars for cohort in cohorts}
    if len(dimensions) != 1:
        raise ValueError(f"Latent dimensions do not match across cohort entries: {sorted(dimensions)}")
    return ad.concat(cohorts, axis=0, join="inner", merge="same", index_unique=None)


def _faiss_neighbors(data: np.ndarray, n_neighbors: int, nlist: int, nprobe: int, gpu: bool) -> tuple[np.ndarray, np.ndarray]:
    try:
        import faiss
    except ImportError as exc:
        raise RuntimeError("The faiss engine requires the faiss Python package") from exc
    quantizer = faiss.IndexFlatL2(data.shape[1])
    index = faiss.IndexIVFFlat(quantizer, data.shape[1], nlist, faiss.METRIC_L2)
    gpu_resources = None
    if gpu:
        if not hasattr(faiss, "StandardGpuResources"):
            raise RuntimeError("GPU backend requires a FAISS build with GPU support")
        # FAISS GPU indexes do not own their resource provider. Keep it alive
        # until all index operations have completed.
        gpu_resources = faiss.StandardGpuResources()
        index = faiss.index_cpu_to_gpu(gpu_resources, 0, index)
    index.train(data)
    index.add(data)
    index.nprobe = nprobe
    distances, indices = index.search(data, n_neighbors + 1)
    return indices[:, 1:].astype(np.int64), distances[:, 1:].astype(np.float32)


def _cpu_umap(data: np.ndarray, indices: np.ndarray, distances: np.ndarray, params: Mapping[str, Any]) -> np.ndarray:
    try:
        import umap
    except ImportError as exc:
        raise RuntimeError("CPU FAISS clustering requires umap-learn") from exc
    model = umap.UMAP(
        n_neighbors=indices.shape[1],
        n_components=2,
        min_dist=float(params["umap_min_dist"]),
        spread=float(params["umap_spread"]),
        random_state=int(params["seed"]),
        precomputed_knn=(indices, distances, None),
    )
    return np.asarray(model.fit_transform(data), dtype=np.float32)


def _cpu_leiden(graph: sparse.csr_matrix, resolution: float, seed: int) -> np.ndarray:
    try:
        import igraph as ig
        import leidenalg
    except ImportError as exc:
        raise RuntimeError("CPU FAISS clustering requires igraph and leidenalg") from exc
    upper = sparse.triu(graph, k=1).tocoo()
    igraph = ig.Graph(n=graph.shape[0], edges=list(zip(upper.row.tolist(), upper.col.tolist())), directed=False)
    igraph.es["weight"] = upper.data.astype(float).tolist()
    partition = leidenalg.find_partition(
        igraph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=float(resolution),
        seed=int(seed),
    )
    return np.asarray(partition.membership).astype(str)


def _gpu_cluster(
    data: np.ndarray,
    indices: np.ndarray,
    distances: np.ndarray,
    graph: sparse.csr_matrix,
    params: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    try:
        import cudf
        import cugraph
        import cupy as cp
        from cuml.manifold import UMAP
        from cupyx.scipy.sparse import coo_matrix
    except ImportError as exc:
        raise RuntimeError("GPU FAISS clustering requires cupy, cudf, cugraph, and cuml") from exc

    upper = sparse.triu(graph, k=1).tocoo()
    edges = cudf.DataFrame(
        {
            "src": upper.row.astype(np.int32),
            "dst": upper.col.astype(np.int32),
            "weight": upper.data.astype(np.float32),
        }
    )
    gpu_graph = cugraph.Graph()
    gpu_graph.from_cudf_edgelist(edges, source="src", destination="dst", edge_attr="weight", renumber=False)
    partitions, _ = cugraph.leiden(
        gpu_graph,
        resolution=float(params["leiden_resolution"]),
        max_iter=1000,
    )
    labels = (
        partitions.to_pandas()
        .set_index("vertex")
        .reindex(np.arange(data.shape[0]))["partition"]
        .astype(str)
        .to_numpy()
    )

    rows = cp.repeat(cp.arange(data.shape[0]), indices.shape[1])
    knn_graph = coo_matrix(
        (
            cp.asarray(distances.ravel()),
            (rows, cp.asarray(indices.ravel())),
        ),
        shape=(data.shape[0], data.shape[0]),
    ).tocsr()
    model = UMAP(
        n_neighbors=indices.shape[1],
        n_components=2,
        min_dist=float(params["umap_min_dist"]),
        spread=float(params["umap_spread"]),
        random_state=int(params["seed"]),
    )
    embedding = model.fit_transform(X=cudf.DataFrame(data), knn_graph=knn_graph)
    if hasattr(embedding, "get"):
        embedding = embedding.get()
    elif hasattr(embedding, "to_numpy"):
        embedding = embedding.to_numpy()
    return labels, np.asarray(embedding)


def cluster_stage(adata: ad.AnnData, cfg: Mapping[str, Any]) -> tuple[ad.AnnData, dict[str, Any]]:
    engine = str(cfg.get("subcluster_engine", "faiss")).lower()
    backend = str(cfg.get("subcluster_backend", "cpu")).lower()
    gpu_postprocess = str(cfg.get("subcluster_gpu_postprocess", "cpu")).lower()
    if engine not in {"faiss", "scanpy"}:
        raise ValueError("subcluster_engine must be faiss or scanpy")
    if backend not in {"cpu", "gpu"}:
        raise ValueError("subcluster_backend must be cpu or gpu")
    if gpu_postprocess not in {"cpu", "rapids"}:
        raise ValueError("subcluster_gpu_postprocess must be cpu or rapids")
    if engine == "scanpy" and backend != "cpu":
        raise ValueError("The scanpy engine only supports the cpu backend")
    if backend == "gpu":
        try:
            import cupy as cp
        except ImportError as exc:
            raise RuntimeError("GPU backend requires cupy and a visible CUDA device") from exc
        try:
            device_count = cp.cuda.runtime.getDeviceCount()
        except Exception as exc:
            raise RuntimeError(
                "GPU backend requested, but CUDA device discovery failed. "
                "Check the Slurm GPU allocation, Singularity --nv, and nvidia-smi on the compute node."
            ) from exc
        if device_count < 1:
            raise RuntimeError("GPU backend requested, but no CUDA device is visible")

    params = {
        "engine": engine,
        "backend": backend,
        "gpu_postprocess": gpu_postprocess if backend == "gpu" else None,
        "seed": int(cfg.get("subcluster_seed", 42)),
        "n_neighbors": int(cfg.get("subcluster_n_neighbors", 30)),
        "leiden_resolution": float(cfg.get("subcluster_leiden_resolution", 0.5)),
        "umap_min_dist": float(cfg.get("subcluster_umap_min_dist", 0.03)),
        "umap_spread": float(cfg.get("subcluster_umap_spread", 0.6)),
    }
    result = adata.copy()
    data = _dense_float32(result.X)
    if engine == "scanpy":
        try:
            import scanpy as sc
        except ImportError as exc:
            raise RuntimeError("The scanpy engine requires scanpy") from exc
        params["n_neighbors"] = min(params["n_neighbors"], result.n_obs - 1)
        sc.pp.neighbors(result, n_neighbors=params["n_neighbors"], use_rep="X", random_state=params["seed"])
        sc.tl.umap(
            result,
            min_dist=params["umap_min_dist"],
            spread=params["umap_spread"],
            random_state=params["seed"],
        )
        sc.tl.leiden(
            result,
            key_added="leiden",
            resolution=params["leiden_resolution"],
            random_state=params["seed"],
            flavor="igraph",
            directed=False,
        )
        return result, params

    n_neighbors, nlist, nprobe = effective_faiss_params(
        result.n_obs,
        params["n_neighbors"],
        cfg.get("subcluster_faiss_nlist"),
        cfg.get("subcluster_faiss_nprobe"),
    )
    params.update({"n_neighbors": n_neighbors, "faiss_nlist": nlist, "faiss_nprobe": nprobe})
    indices, distances = _faiss_neighbors(data, n_neighbors, nlist, nprobe, backend == "gpu")
    graph = symmetric_affinity_graph(indices, distances, result.n_obs)
    if backend == "gpu" and gpu_postprocess == "rapids":
        labels, embedding = _gpu_cluster(data, indices, distances, graph, params)
    else:
        labels = _cpu_leiden(graph, params["leiden_resolution"], params["seed"])
        embedding = _cpu_umap(data, indices, distances, params)
    result.obs["leiden"] = pd.Categorical(labels)
    result.obsm["X_umap"] = np.asarray(embedding, dtype=np.float32)
    result.obsp["connectivities"] = graph
    return result, params


def write_diagnostics(stage_name: str, result: ad.AnnData, params: Mapping[str, Any], output_dir: pathlib.Path) -> None:
    write_compatible_h5ad(result, output_dir / f"query_latent_tumor_{stage_name}.h5ad")
    pd.DataFrame(
        {"index": result.obs["source_index"].astype(str), "leiden": result.obs["leiden"].astype(str)}
    ).to_csv(output_dir / f"atlas_{stage_name}_leiden.csv", index=False)
    pd.DataFrame(
        result.obsm["X_umap"],
        index=result.obs["source_index"].astype(str),
        columns=["UMAP1", "UMAP2"],
    ).to_csv(output_dir / f"atlas_{stage_name}_umap.csv")
    write_json(output_dir / f"atlas_{stage_name}_uparams.json", params)

    counts = pd.crosstab(result.obs["leiden"].astype(str), result.obs["cell_type"].astype(str))
    percentages = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0) * 100
    counts.to_csv(output_dir / f"contingency_{stage_name}_counts.csv")
    percentages.to_csv(output_dir / f"contingency_{stage_name}_percent.csv")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    embedding = result.obsm["X_umap"]
    for axis, key in zip(axes, ["leiden", "dataset"]):
        labels = result.obs[key].astype(str)
        for label in sorted(labels.unique()):
            mask = labels == label
            axis.scatter(embedding[mask, 0], embedding[mask, 1], s=2, alpha=0.7, label=label)
        axis.set_title(f"{stage_name}: {key}")
        axis.set_xlabel("UMAP1")
        axis.set_ylabel("UMAP2")
        axis.legend(markerscale=4, fontsize=6)
    fig.tight_layout()
    fig.savefig(output_dir / f"umap_{stage_name}.png", dpi=180)
    plt.close(fig)

    fig_width = max(10, min(30, 0.45 * max(1, percentages.shape[1])))
    fig_height = max(5, min(24, 0.45 * max(1, percentages.shape[0])))
    fig, axis = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(percentages, cmap="Blues", ax=axis)
    axis.set_title(f"{stage_name}: percent of Leiden cluster")
    fig.tight_layout()
    fig.savefig(output_dir / f"contingency_{stage_name}.png", dpi=180)
    plt.close(fig)


def run(cfg: Mapping[str, Any]) -> pathlib.Path:
    output_dir = pathlib.Path(cfg.get("subcluster_output_dir", "subcluster_results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    merged = load_cohort(cfg)
    merged.obs_names_make_unique()
    write_compatible_h5ad(merged, output_dir / "query_latent.h5ad")

    selected = select_subcluster_cells(
        merged.obs,
        tumor_pattern=str(cfg.get("subcluster_tumor_pattern", "Tumor")),
        epithelial_types=cfg.get("subcluster_epithelial_types") or DEFAULT_EPITHELIAL_TYPES,
    )
    stage_masks = split_stage_masks(merged.obs)
    qc: dict[str, Any] = {
        "n_obs_merged": merged.n_obs,
        "n_latent_dimensions": merged.n_vars,
        "n_selected_tumor_epithelial": int(selected.sum()),
        "n_unknown_stage": int(stage_masks["unknown"].sum()),
        "unknown_stage_by_dataset": (
            merged.obs.loc[stage_masks["unknown"]]
            .groupby("dataset", observed=True)
            .size()
            .astype(int)
            .to_dict()
        ),
        "stages": {},
    }
    for stage_name in ("early", "late"):
        mask = selected & stage_masks[stage_name]
        if int(mask.sum()) < 3:
            raise ValueError(f"{stage_name} subclustering selected only {int(mask.sum())} cells")
        stage_result, params = cluster_stage(merged[mask].copy(), cfg)
        write_diagnostics(stage_name, stage_result, params, output_dir)
        qc["stages"][stage_name] = {
            "n_obs": stage_result.n_obs,
            "n_clusters": int(stage_result.obs["leiden"].nunique()),
            "datasets": stage_result.obs.groupby("dataset", observed=True).size().astype(int).to_dict(),
        }
    write_json(output_dir / "qc_summary.json", qc)
    write_json(
        output_dir / "versions.json",
        {
            name: package_version(name)
            for name in [
                "anndata",
                "scanpy",
                "faiss-cpu",
                "faiss-gpu",
                "igraph",
                "leidenalg",
                "umap-learn",
                "cudf",
                "cugraph",
                "cuml",
                "cupy",
            ]
        },
    )
    return output_dir


def apply_mappings(raw_dir: str | pathlib.Path, cfg: Mapping[str, Any], output_dir: str | pathlib.Path) -> pathlib.Path:
    raw_dir = pathlib.Path(raw_dir)
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    merged = ad.read_h5ad(raw_dir / "query_latent.h5ad")
    merged.obs["cell_type_adjusted"] = merged.obs["cell_type"].astype(str)
    status: dict[str, Any] = {"applied": [], "skipped": []}
    for stage_name, config_key in [("early", "early_mapping_json"), ("late", "late_mapping_json")]:
        mapping_path = cfg.get(config_key)
        if not mapping_path:
            status["skipped"].append(stage_name)
            continue
        stage_result = ad.read_h5ad(raw_dir / f"query_latent_tumor_{stage_name}.h5ad")
        mapping = read_structured(mapping_path)
        mapping = mapping.get("mapping", mapping)
        validate_mapping(stage_result.obs["leiden"].astype(str), mapping)
        write_json(output_dir / f"mapping_{stage_name}_leiden.json", mapping)
        apply_mapping_exports(
            merged,
            stage_name,
            stage_result.obs_names,
            stage_result.obs["leiden"].astype(str),
            mapping,
            output_dir,
        )
        mapped = pd.Series(
            [str(mapping[str(label)]) for label in stage_result.obs["leiden"].astype(str)],
            index=stage_result.obs_names,
            dtype="object",
        )
        merged.obs.loc[mapped.index, "cell_type_adjusted"] = mapped
        status["applied"].append(stage_name)
    write_compatible_h5ad(merged, output_dir / "query_latent_adjusted.h5ad")
    write_json(output_dir / "mapping_status.json", status)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare and subcluster the LUCA latent cohort")
    parser.add_argument("--config", required=True, help="YAML or JSON subclustering configuration")
    parser.add_argument("--apply-mappings-from", help="Raw subcluster results directory")
    parser.add_argument("--mapping-output-dir", default="mapping_results")
    args = parser.parse_args()
    cfg = read_structured(args.config)
    if args.apply_mappings_from:
        apply_mappings(args.apply_mappings_from, cfg, args.mapping_output_dir)
    else:
        run(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
