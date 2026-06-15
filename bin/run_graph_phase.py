#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import subprocess
import sys
from collections.abc import Mapping
from typing import Any

import anndata as ad
import igraph as ig
import leidenalg
import matplotlib
import networkx as nx
import numpy as np
import pandas as pd
import yaml
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
from sklearn.metrics import silhouette_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EARLY = {"I", "II"}
LATE = {"III", "IV", "III or IV"}


def read_structured(path: str | pathlib.Path) -> dict[str, Any]:
    path = pathlib.Path(path)
    with path.open() as handle:
        return yaml.safe_load(handle) if path.suffix.lower() in {".yaml", ".yml"} else json.load(handle)


def write_json(path: pathlib.Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str))


def canonical_stage(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().upper().replace("_", " ").replace("-", " ")
    if text in {"III OR IV", "III/IV", "III IV"}:
        return "III or IV"
    for stage in ("IV", "III", "II", "I"):
        if text.startswith(stage):
            return stage
    return None


def load_cells(cfg: Mapping[str, Any]) -> pd.DataFrame:
    input_path = pathlib.Path(cfg["graph_input"])
    if input_path.suffix.lower() == ".csv":
        obs = pd.read_csv(input_path)
    else:
        obs = ad.read_h5ad(input_path, backed="r").obs.copy()
    label_key = cfg.get("graph_cell_type_key", "cell_type_adjusted")
    if label_key not in obs and label_key == "cell_type_adjusted" and "cell_type" in obs:
        label_key = "cell_type"
    keys = {
        "cell_type": label_key,
        "sample": cfg.get("graph_sample_key", "batch"),
        "dataset": cfg.get("graph_dataset_key", "dataset"),
        "stage": cfg.get("graph_stage_key", "stage"),
    }
    missing = [source for source in keys.values() if source not in obs]
    if missing:
        raise ValueError(f"Graph input {input_path} is missing columns: {missing}")
    result = obs[[keys["cell_type"], keys["sample"], keys["dataset"], keys["stage"]]].copy()
    result.columns = ["cell_type_adjusted", "sample", "dataset", "stage"]
    result = result.dropna(subset=["cell_type_adjusted", "sample"])
    result["stage"] = result["stage"].map(canonical_stage)
    for column in ["cell_type_adjusted", "sample", "dataset"]:
        result[column] = result[column].astype(str)
    return result


def prepare_global(cfg: Mapping[str, Any], output: pathlib.Path) -> None:
    cells = load_cells(cfg)
    global_dir = output / "global"
    global_dir.mkdir(parents=True, exist_ok=True)
    unknown = cells[cells["stage"].isna()].groupby("dataset", observed=True).size()
    summary: dict[str, Any] = {"input_cells": len(cells), "unknown_stage_cells": int(unknown.sum()), "stages": {}}
    min_count = int(cfg.get("graph_min_cells_per_type_sample", 2))
    for name, stages in [("early", EARLY), ("late", LATE)]:
        selected = cells[cells["stage"].isin(stages)]
        groups = (
            selected.groupby(["cell_type_adjusted", "sample", "dataset"], observed=True)
            .size().rename("count").reset_index()
        )
        groups.to_csv(global_dir / f"groups_{name}.csv", index=False)
        filtered = groups[groups["count"] >= min_count]
        matrix = filtered.pivot_table(
            index="cell_type_adjusted", columns="sample", values="count", aggfunc="sum", fill_value=0
        )
        matrix.index.name = "types"
        matrix.to_csv(global_dir / f"matrix_{name}.txt", sep="\t")
        sample_meta = groups[["sample", "dataset"]].drop_duplicates("sample")
        sample_meta.to_csv(global_dir / f"samples_{name}.csv", index=False)
        summary["stages"][name] = {
            "cells": len(selected), "samples": matrix.shape[1], "cell_types": matrix.shape[0]
        }
    write_json(global_dir / "qc_summary.json", summary)


def cluster_samples(matrix: pd.DataFrame, k: int, resolution: float, seed: int, target_ecotypes: int | None = None, optimize_resolution: bool = False) -> tuple[pd.Series, nx.Graph]:
    abundance = matrix.T.astype(float)
    abundance = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    n = len(abundance)
    if n < 2:
        raise ValueError(f"Ecotype clustering requires at least two samples; got {n}")
    k = min(max(1, k), n - 1)
    distances = cdist(abundance, abundance)
    graph = nx.Graph()
    graph.add_nodes_from(abundance.index.astype(str))
    for i, sample in enumerate(abundance.index.astype(str)):
        for j in np.argsort(distances[i])[1 : k + 1]:
            graph.add_edge(sample, str(abundance.index[j]))
    igraph = ig.Graph.TupleList(graph.edges(), directed=False)
    missing = sorted(set(graph.nodes) - set(igraph.vs["name"]))
    igraph.add_vertices(missing)
    
    best_partition = None
    if target_ecotypes is not None and target_ecotypes > 0:
        low, high = 0.01, 5.0
        best_diff = float('inf')
        for _ in range(20):
            mid = (low + high) / 2
            partition = leidenalg.find_partition(
                igraph, leidenalg.RBConfigurationVertexPartition,
                resolution_parameter=mid, seed=seed,
            )
            n_clusters = len(set(partition.membership))
            if n_clusters == target_ecotypes:
                best_partition = partition
                break
            elif n_clusters < target_ecotypes:
                low = mid
            else:
                high = mid
            if abs(n_clusters - target_ecotypes) < best_diff:
                best_diff = abs(n_clusters - target_ecotypes)
                best_partition = partition
        if best_partition is None:
            best_partition = leidenalg.find_partition(
                igraph, leidenalg.RBConfigurationVertexPartition,
                resolution_parameter=resolution, seed=seed,
            )
    elif optimize_resolution:
        best_score = -1
        resolutions = np.linspace(0.05, 2.0, 40)
        for res in resolutions:
            partition = leidenalg.find_partition(
                igraph, leidenalg.RBConfigurationVertexPartition,
                resolution_parameter=res, seed=seed,
            )
            labels = np.array(partition.membership)
            n_clusters = len(set(labels))
            if 1 < n_clusters < len(abundance):
                try:
                    score = silhouette_score(distances, labels, metric='precomputed')
                    if score > best_score:
                        best_score = score
                        best_partition = partition
                except ValueError:
                    pass
        if best_partition is None:
            best_partition = leidenalg.find_partition(
                igraph, leidenalg.RBConfigurationVertexPartition,
                resolution_parameter=resolution, seed=seed,
            )
    else:
        best_partition = leidenalg.find_partition(
            igraph, leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution, seed=seed,
        )

    membership = pd.Series(best_partition.membership, index=igraph.vs["name"], name="membership")
    return membership.reindex(abundance.index.astype(str)).astype(int), graph


def run_ecotype(cfg: Mapping[str, Any], output: pathlib.Path) -> None:
    global_dir = pathlib.Path(cfg.get("graph_global_dir") or output / "global")
    ecotype_dir = output / "ecotype"
    ecotype_dir.mkdir(parents=True, exist_ok=True)
    for stage in ("early", "late"):
        matrix = pd.read_csv(global_dir / f"matrix_{stage}.txt", sep="\t", index_col=0)
        metadata = pd.read_csv(global_dir / f"samples_{stage}.csv").set_index("sample")
        target_ecotypes = cfg.get("ecotype_target_k")
        target_ecotypes = int(target_ecotypes) if target_ecotypes is not None else None
        optimize_resolution = bool(cfg.get("ecotype_optimize_resolution", False))
        membership, graph = cluster_samples(
            matrix, int(cfg.get("ecotype_k", 4)), float(cfg.get("ecotype_resolution", 0.2)),
            int(cfg.get("graph_seed", 42)), target_ecotypes=target_ecotypes,
            optimize_resolution=optimize_resolution
        )
        result = membership.rename_axis("sample").reset_index()
        result.insert(1, "dataset", result["sample"].map(metadata["dataset"]).fillna("unknown"))
        result.to_csv(ecotype_dir / f"membership_{stage}.csv", index=False)
        abundance = matrix.T.astype(float)
        abundance = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        abundance.insert(0, "membership", membership.reindex(abundance.index).astype(int))
        abundance.insert(0, "dataset", abundance.index.map(metadata["dataset"]).fillna("unknown"))
        abundance.to_csv(ecotype_dir / f"abundance_{stage}.csv")
        nx.write_graphml(graph, ecotype_dir / f"sample_knn_{stage}.graphml")
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(graph, seed=int(cfg.get("graph_seed", 42)))
        nx.draw_networkx(graph, pos, node_size=80, font_size=5, node_color=membership.reindex(graph.nodes).values)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(ecotype_dir / f"sample_clusters_{stage}.png", dpi=160)
        plt.close()


def write_aracne_inputs(matrix: pd.DataFrame, directory: pathlib.Path, prefix: str) -> tuple[pathlib.Path, pathlib.Path]:
    directory.mkdir(parents=True, exist_ok=True)
    matrix = matrix.loc[matrix.sum(axis=1) > 0].copy()
    matrix.index.name = "types"
    matrix_path = directory / f"matrix_{prefix}.txt"
    tf_path = directory / f"tfs_{prefix}.txt"
    matrix.to_csv(matrix_path, sep="\t")
    tf_path.write_text("\n".join(matrix.index.astype(str)) + "\n")
    return matrix_path, tf_path


def run_aracne(matrix: pd.DataFrame, directory: pathlib.Path, prefix: str, cfg: Mapping[str, Any]) -> pathlib.Path | None:
    matrix_path, tf_path = write_aracne_inputs(matrix, directory, prefix)
    if bool(cfg.get("aracne_dry_run", False)):
        return None
    jar = pathlib.Path(cfg["aracne_jar"])
    if not jar.exists():
        raise FileNotFoundError(f"ARACNe jar does not exist: {jar}")
    java = str(cfg.get("aracne_java", "java"))
    if not pathlib.Path(java).exists() and not shutil.which(java):
        raise RuntimeError(f"Java executable is required to run ARACNe; not found: {java}")
    base = [java, f"-Xmx{cfg.get('aracne_java_memory', '8G')}", "-jar", str(jar.resolve()),
            "-e", str(matrix_path.resolve()), "-o", str(directory.resolve()), "--tfs", str(tf_path.resolve()),
            "--pvalue", str(cfg.get("aracne_pvalue", "1E-8"))]
    subprocess.run(base + ["--seed", "1", "--calculateThreshold"], check=True)
    for seed in range(1, int(cfg.get("aracne_bootstraps", 500)) + 1):
        subprocess.run(base + ["--seed", str(seed)], check=True)
    subprocess.run([java, f"-Xmx{cfg.get('aracne_java_memory', '8G')}", "-jar", str(jar.resolve()),
                    "-o", str(directory.resolve()), "--consolidate"], check=True)
    network = directory / "network.txt"
    destination = directory / f"net_{prefix}.txt"
    lines = network.read_text().splitlines()
    destination.write_text("\n".join(lines[1:]) + ("\n" if len(lines) > 1 else ""))
    return destination


def annotate_network(network_path: pathlib.Path, matrix: pd.DataFrame, output: pathlib.Path) -> None:
    try:
        edges = pd.read_csv(network_path, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        output.write_text("")
        return
    edges.columns = ["source", "target", "MI", "p", *[f"extra_{i}" for i in range(max(0, edges.shape[1] - 4))]]
    abundance = matrix.T.astype(float)
    abundance = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    correlations, pvalues = [], []
    for row in edges.itertuples(index=False):
        corr, pvalue = pearsonr(abundance[str(row.source)], abundance[str(row.target)])
        correlations.append(0.0 if np.isnan(corr) else corr)
        pvalues.append(1.0 if np.isnan(pvalue) else pvalue)
    edges["pearson"] = correlations
    edges["pearson_pvalue"] = pvalues
    edges.to_csv(output, sep="\t", index=False)


def run_global_aracne(cfg: Mapping[str, Any], output: pathlib.Path) -> None:
    global_dir = pathlib.Path(cfg.get("graph_global_dir") or output / "global")
    for stage in ("early", "late"):
        matrix = pd.read_csv(global_dir / f"matrix_{stage}.txt", sep="\t", index_col=0)
        network = run_aracne(matrix, global_dir / "aracne" / stage, stage, cfg)
        if network:
            annotate_network(network, matrix, network.with_name(f"net_{stage}_MI_pearson.txt"))


def run_ecotype_aracne(cfg: Mapping[str, Any], output: pathlib.Path) -> None:
    global_dir = pathlib.Path(cfg.get("graph_global_dir") or output / "global")
    ecotype_dir = pathlib.Path(cfg.get("ecotype_dir") or output / "ecotype")
    cluster_root = output / "ecotype_graphs"
    min_samples = int(cfg.get("aracne_min_samples", 5))
    min_types = int(cfg.get("aracne_min_nonzero_types", 3))
    for stage in ("early", "late"):
        matrix = pd.read_csv(global_dir / f"matrix_{stage}.txt", sep="\t", index_col=0)
        membership = pd.read_csv(ecotype_dir / f"membership_{stage}.csv")
        for cluster, rows in membership.groupby("membership", observed=True):
            samples = [sample for sample in rows["sample"].astype(str) if sample in matrix.columns]
            subset = matrix.loc[:, samples]
            subset = subset.loc[subset.sum(axis=1) > 0]
            directory = cluster_root / stage / f"cluster_{int(cluster):02d}"
            directory.mkdir(parents=True, exist_ok=True)
            write_json(directory / "qc.json", {"samples": len(samples), "cell_types": len(subset)})
            if len(samples) < min_samples or len(subset) < min_types:
                continue
            prefix = f"{stage}_cluster{int(cluster)}"
            network = run_aracne(subset, directory, prefix, cfg)
            if network:
                annotate_network(network, subset, directory / f"net_{prefix}_MI_pearson.txt")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--step", choices=["prepare", "global-aracne", "ecotype", "ecotype-aracne", "all"], default="all")
    args = parser.parse_args()
    cfg = read_structured(args.config)
    output = pathlib.Path(cfg.get("graph_output_dir", "graph_results"))
    output.mkdir(parents=True, exist_ok=True)
    steps = ["prepare", "global-aracne", "ecotype", "ecotype-aracne"] if args.step == "all" else [args.step]
    for step in steps:
        {"prepare": prepare_global, "global-aracne": run_global_aracne, "ecotype": run_ecotype,
         "ecotype-aracne": run_ecotype_aracne}[step](cfg, output)
    write_json(output / "run_config.json", cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
