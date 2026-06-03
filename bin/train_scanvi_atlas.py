#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import json
import os
import pathlib
import random
import shutil
import sys
import time
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import yaml
from lightning.pytorch.loggers import TensorBoardLogger


def read_config(path: pathlib.Path) -> dict[str, Any]:
    with path.open() as handle:
        if path.suffix.lower() in {".yaml", ".yml"}:
            return yaml.safe_load(handle)
        return json.load(handle)


def as_jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, pathlib.Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(k): as_jsonable(v) for k, v in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [as_jsonable(v) for v in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return str(value)


def write_json(path: pathlib.Path, data: Any) -> None:
    path.write_text(json.dumps(as_jsonable(data), indent=2, sort_keys=True))


def package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def filter_obs_values(adata: sc.AnnData, key: str | None, values: Any, keep: bool) -> sc.AnnData:
    values = as_list(values)
    if not key or not values:
        return adata
    if key not in adata.obs:
        raise KeyError(f"adata.obs has no column {key!r}")
    mask = adata.obs[key].isin(values)
    return adata[mask if keep else ~mask].copy()


def maybe_split_batch_by_sample(adata: sc.AnnData, cfg: Mapping[str, Any]) -> None:
    dataset = cfg.get("split_batch_dataset")
    batch_key = cfg["batch_key"]
    sample_key = cfg.get("sample_key", "sample")
    if not dataset:
        return
    if batch_key not in adata.obs:
        raise KeyError(f"adata.obs has no batch_key column {batch_key!r}")
    if sample_key not in adata.obs:
        raise KeyError(f"adata.obs has no sample_key column {sample_key!r}")
    batch = adata.obs[batch_key].astype(str)
    sample = adata.obs[sample_key].astype(str)
    adata.obs[batch_key] = np.where(batch == str(dataset), batch + "_" + sample, batch)


def prepare_counts(adata: sc.AnnData, cfg: Mapping[str, Any]) -> sc.AnnData:
    counts_layer = cfg.get("counts_layer")
    if counts_layer:
        if counts_layer not in adata.layers:
            raise KeyError(f"adata.layers has no counts_layer {counts_layer!r}")
        adata.X = adata.layers[counts_layer].copy()
    if cfg.get("round_counts", False):
        if hasattr(adata.X, "data"):
            adata.X.data = np.round(adata.X.data)
        else:
            adata.X = np.round(adata.X)
    for layer in as_list(cfg.get("drop_layers")):
        if layer in adata.layers:
            del adata.layers[layer]
    if cfg.get("drop_raw", True):
        adata.raw = None
    return adata


def subset_hvg(adata: sc.AnnData, cfg: Mapping[str, Any]) -> sc.AnnData:
    mode = cfg.get("hvg_mode", "compute")
    if mode in {None, "none", "all"}:
        return adata
    if mode == "precomputed":
        column = cfg.get("hvg_column", "is_highly_variable")
        if column not in adata.var:
            raise KeyError(f"adata.var has no hvg_column {column!r}")
        mask = adata.var[column].astype(str).str.lower().isin(["true", "1", "yes"])
        if not bool(mask.sum()):
            raise ValueError(f"hvg_column {column!r} did not select any genes")
        return adata[:, mask.to_numpy()].copy()
    if mode != "compute":
        raise ValueError(f"Unsupported hvg_mode {mode!r}; use compute, precomputed, or none")

    batch_key = cfg.get("hvg_batch_key") or cfg.get("batch_key")
    kwargs = {
        "flavor": cfg.get("hvg_flavor", "seurat_v3"),
        "n_top_genes": int(cfg.get("hvg_n_top_genes", 6000)),
        "inplace": False,
    }
    if batch_key:
        kwargs["batch_key"] = batch_key
    hvg = sc.pp.highly_variable_genes(adata, **kwargs)
    adata.var["is_highly_variable"] = hvg["highly_variable"].to_numpy()
    return adata[:, hvg["highly_variable"].to_numpy()].copy()


def load_and_prepare_adata(cfg: Mapping[str, Any]) -> sc.AnnData:
    adata = sc.read_h5ad(cfg["input_h5ad"])
    adata = filter_obs_values(adata, cfg.get("origin_key"), cfg.get("origin_values"), keep=True)
    adata = filter_obs_values(adata, cfg.get("stage_key"), cfg.get("stages"), keep=True)
    adata = filter_obs_values(adata, cfg.get("study_key"), cfg.get("exclude_studies"), keep=False)
    adata = prepare_counts(adata.copy(), cfg)
    maybe_split_batch_by_sample(adata, cfg)
    max_cells = cfg.get("max_cells")
    if max_cells and adata.n_obs > int(max_cells):
        sc.pp.subsample(adata, n_obs=int(max_cells), random_state=int(cfg.get("seed", 0)), copy=False)
    adata = subset_hvg(adata, cfg)
    max_genes = cfg.get("max_genes")
    if max_genes and adata.n_vars > int(max_genes):
        adata = adata[:, : int(max_genes)].copy()
    return adata


def history_to_frame(history: Any) -> pd.DataFrame:
    if history is None:
        return pd.DataFrame()
    if isinstance(history, pd.DataFrame):
        return history.reset_index(drop=False)
    if not isinstance(history, Mapping):
        return pd.DataFrame()
    columns: dict[str, pd.Series] = {}
    for key, value in history.items():
        if isinstance(value, pd.DataFrame):
            for column in value.columns:
                columns[f"{key}.{column}"] = pd.Series(value[column].to_numpy())
        elif isinstance(value, pd.Series):
            columns[str(key)] = value.reset_index(drop=True)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
            columns[str(key)] = pd.Series(value)
    return pd.DataFrame(columns)


def write_history(path: pathlib.Path, history: Any) -> None:
    frame = history_to_frame(history)
    if not frame.empty:
        frame.to_csv(path, index=False)


def coerce_numeric_strings(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: coerce_numeric_strings(item) for key, item in value.items()}
    if isinstance(value, list):
        return [coerce_numeric_strings(item) for item in value]
    if isinstance(value, str):
        try:
            return float(value) if any(char in value.lower() for char in [".", "e"]) else int(value)
        except ValueError:
            return value
    return value


def tensorboard_logger(cfg: Mapping[str, Any], run_dir: pathlib.Path, stage: str) -> TensorBoardLogger | bool:
    if not cfg.get("tensorboard_logging", True):
        return False
    log_dir = pathlib.Path(cfg.get("tensorboard_log_dir") or (run_dir / "tensorboard"))
    log_dir.mkdir(parents=True, exist_ok=True)
    return TensorBoardLogger(
        save_dir=str(log_dir),
        name=stage,
        version="version_0",
        default_hp_metric=False,
    )


def build_latent(scanvi_model: scvi.model.SCANVI, adata: sc.AnnData, cfg: Mapping[str, Any]) -> sc.AnnData:
    labels_key = cfg["labels_key"]
    batch_key = cfg["batch_key"]
    latent = sc.AnnData(scanvi_model.get_latent_representation())
    latent.obs_names = adata.obs_names.copy()
    latent.obs["cell_type"] = adata.obs[labels_key].astype(str).to_numpy()
    latent.obs["batch"] = adata.obs[batch_key].astype(str).to_numpy()
    latent.obs["predictions"] = scanvi_model.predict()
    latent.obs["is_correct"] = latent.obs["predictions"].astype(str) == latent.obs["cell_type"].astype(str)
    for column in as_list(cfg.get("latent_obs_columns")):
        if column in adata.obs and column not in latent.obs:
            latent.obs[column] = adata.obs[column].to_numpy()
    if cfg.get("collapse_wu_dataset", True):
        latent.obs["dataset_collapsed"] = latent.obs["batch"].map(
            lambda value: "Wu_Zhou_2021" if "Wu_Zhou_2021" in str(value) else str(value)
        )
    return latent


def fix_arrow_string_arrays(adata: sc.AnnData) -> None:
    """Convert pandas ArrowStringArray columns to regular Python strings for h5py compatibility."""
    for col in adata.obs.columns:
        try:
            # Check if column is ArrowStringArray
            if hasattr(adata.obs[col], 'dtype') and 'string' in str(adata.obs[col].dtype):
                # Convert to regular Python string dtype
                adata.obs[col] = adata.obs[col].astype("object").astype(str)
        except Exception as exc:
            print(f"Warning: could not convert column {col}: {exc}", file=sys.stderr)


def maybe_compute_neighbors(latent: sc.AnnData, cfg: Mapping[str, Any]) -> None:
    if not cfg.get("compute_neighbors", True):
        return
    numba_threads = cfg.get("numba_threads")
    if numba_threads:
        try:
            import numba

            numba.set_num_threads(int(numba_threads))
        except Exception as exc:
            print(f"Warning: could not set numba_threads={numba_threads}: {exc}", file=sys.stderr)
    sc.pp.neighbors(latent, n_neighbors=int(cfg.get("n_neighbors", 8)))
    if cfg.get("compute_leiden", True):
        sc.tl.leiden(latent)
    if cfg.get("compute_umap", True):
        sc.tl.umap(latent)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=pathlib.Path)
    args = parser.parse_args()

    cfg = read_config(args.config)
    cwd = pathlib.Path.cwd()
    run_dir = pathlib.Path(cfg["run_dir"])
    run_dir.mkdir(parents=True, exist_ok=True)

    seed = int(cfg.get("seed", 0))
    random.seed(seed)
    np.random.seed(seed)
    scvi.settings.seed = seed
    scvi.settings.dl_num_workers = int(cfg.get("dl_num_workers", 8))
    if cfg.get("torch_multiprocessing_sharing_strategy"):
        torch.multiprocessing.set_sharing_strategy(cfg["torch_multiprocessing_sharing_strategy"])
    torch.set_float32_matmul_precision(cfg.get("float32_matmul_precision", "high"))

    versions = {
        "python": sys.version,
        "anndata": package_version("anndata"),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scanpy": package_version("scanpy"),
        "scarches": package_version("scarches"),
        "scvi": scvi.__version__,
        "torch": torch.__version__,
        "torch_cuda_available": torch.cuda.is_available(),
        "torch_cuda_device_count": torch.cuda.device_count(),
        "torch_cuda_device_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
    }
    write_json(cwd / "versions.json", versions)
    write_json(run_dir / "versions.json", versions)

    started = time.time()
    adata = load_and_prepare_adata(cfg)
    if cfg.get("save_prepared_h5ad", False):
        prepared_path = run_dir / "prepared_atlas.h5ad"
        adata.write_h5ad(prepared_path)
        shutil.copy2(prepared_path, cwd / "prepared_atlas.h5ad")

    if cfg["batch_key"] not in adata.obs:
        raise KeyError(f"adata.obs has no batch_key column {cfg['batch_key']!r}")
    if cfg["labels_key"] not in adata.obs:
        raise KeyError(f"adata.obs has no labels_key column {cfg['labels_key']!r}")

    scvi.model.SCVI.setup_anndata(adata, labels_key=cfg["labels_key"], batch_key=cfg["batch_key"])

    scvi_path = run_dir / "scvi_model"
    if cfg.get("reuse_scvi_model_if_exists", True) and scvi_path.exists():
        print(f"Loading existing SCVI model from {scvi_path}")
        scvi_model = scvi.model.SCVI.load(scvi_path, adata=adata)
    else:
        scvi_model = scvi.model.SCVI(adata, **coerce_numeric_strings(cfg.get("scvi_model_params") or {}))
        scvi_model.train(
            early_stopping=bool(cfg.get("scvi_early_stopping", True)),
            **coerce_numeric_strings(cfg.get("scvi_train_params") or {}),
            early_stopping_monitor=(cfg.get("scvi_early_stopping_kwargs") or {}).get("early_stopping_monitor", "elbo_validation"),
            early_stopping_patience=(cfg.get("scvi_early_stopping_kwargs") or {}).get("early_stopping_patience", 10),
            early_stopping_min_delta=(cfg.get("scvi_early_stopping_kwargs") or {}).get("early_stopping_min_delta", 0.1),
            plan_kwargs=coerce_numeric_strings(cfg.get("scvi_plan_kwargs") or None),
            logger=tensorboard_logger(cfg, run_dir, "scvi"),
            **coerce_numeric_strings(cfg.get("scvi_trainer_kwargs") or {}),
        )
        write_history(cwd / "scvi_history.csv", scvi_model.history)
        write_history(run_dir / "scvi_history.csv", scvi_model.history)
        scvi_model.save(scvi_path, overwrite=True)

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=cfg.get("unlabeled_category", "Unknown"),
    )
    scanvi_model.train(
        early_stopping=bool(cfg.get("scanvi_early_stopping", True)),
        **(cfg.get("scanvi_train_params") or {}),
        early_stopping_monitor=(cfg.get("scanvi_early_stopping_kwargs") or {}).get(
            "early_stopping_monitor", "elbo_validation"
        ),
        early_stopping_patience=(cfg.get("scanvi_early_stopping_kwargs") or {}).get("early_stopping_patience", 40),
        early_stopping_min_delta=(cfg.get("scanvi_early_stopping_kwargs") or {}).get("early_stopping_min_delta", 0.7),
        plan_kwargs=coerce_numeric_strings(cfg.get("scanvi_plan_kwargs") or None),
        logger=tensorboard_logger(cfg, run_dir, "scanvi"),
        **coerce_numeric_strings(cfg.get("scanvi_trainer_kwargs") or {}),
    )
    write_history(cwd / "scanvi_history.csv", scanvi_model.history)
    write_history(run_dir / "scanvi_history.csv", scanvi_model.history)

    scanvi_path = run_dir / "scanvi_model"
    scanvi_model.save(scanvi_path, overwrite=True)

    latent = build_latent(scanvi_model, adata, cfg)
    maybe_compute_neighbors(latent, cfg)
    fix_arrow_string_arrays(latent)
    latent.write_h5ad(cwd / "ref_latent.h5ad")
    latent.write_h5ad(run_dir / "ref_latent.h5ad")

    accuracy = float(latent.obs["is_correct"].mean()) if latent.n_obs else float("nan")
    metrics = {
        "experiment_name": cfg["experiment_name"],
        "elapsed_seconds_total": time.time() - started,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "batch_key": cfg["batch_key"],
        "labels_key": cfg["labels_key"],
        "accuracy": accuracy,
        "scvi_model": str(scvi_path),
        "scanvi_model": str(scanvi_path),
        "ref_latent": str(run_dir / "ref_latent.h5ad"),
    }
    write_json(cwd / "metrics.json", metrics)
    write_json(run_dir / "metrics.json", metrics)

    del adata, latent, scvi_model, scanvi_model
    gc.collect()
    print(f"Finished {cfg['experiment_name']}; results in {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
