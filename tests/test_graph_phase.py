import importlib.util
import pathlib

import anndata as ad
import numpy as np
import pandas as pd

MODULE_PATH = pathlib.Path(__file__).parents[1] / "bin" / "run_graph_phase.py"
SPEC = importlib.util.spec_from_file_location("run_graph_phase", MODULE_PATH)
graph = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(graph)


def test_canonical_stage():
    assert graph.canonical_stage("IA1") == "I"
    assert graph.canonical_stage("IIIB") == "III"
    assert graph.canonical_stage("III or IV") == "III or IV"
    assert graph.canonical_stage("unknown") is None


def test_prepare_global_builds_stage_matrices(tmp_path):
    obs = pd.DataFrame({
        "cell_type_adjusted": ["A", "A", "B", "B", "A", "A"],
        "batch": ["s1", "s1", "s1", "s1", "s2", "s2"],
        "dataset": ["D"] * 6,
        "stage": ["I", "I", "I", "I", "III", "III"],
    })
    path = tmp_path / "cells.h5ad"
    ad.AnnData(X=np.ones((6, 1)), obs=obs).write_h5ad(path)
    graph.prepare_global({"graph_input": str(path), "graph_min_cells_per_type_sample": 2}, tmp_path / "out")
    early = pd.read_csv(tmp_path / "out/global/matrix_early.txt", sep="\t", index_col=0)
    late = pd.read_csv(tmp_path / "out/global/matrix_late.txt", sep="\t", index_col=0)
    assert early.loc["A", "s1"] == 2
    assert early.loc["B", "s1"] == 2
    assert late.loc["A", "s2"] == 2


def test_ecotype_clustering_is_bounded_and_aligned():
    matrix = pd.DataFrame(
        {"s1": [10, 0], "s2": [9, 1], "s3": [0, 10], "s4": [1, 9]},
        index=["A", "B"],
    )
    membership, network = graph.cluster_samples(matrix, k=99, resolution=0.2, seed=42)
    assert membership.index.tolist() == ["s1", "s2", "s3", "s4"]
    assert network.number_of_nodes() == 4


def test_annotate_network_adds_pearson_columns(tmp_path):
    matrix = pd.DataFrame({"s1": [1, 4], "s2": [2, 3], "s3": [3, 2], "s4": [4, 1]}, index=["A", "B"])
    network = tmp_path / "net.txt"
    network.write_text("A\tB\t0.5\t0.01\n")
    output = tmp_path / "annotated.txt"
    graph.annotate_network(network, matrix, output)
    annotated = pd.read_csv(output, sep="\t")
    assert annotated.loc[0, "pearson"] < -0.9
    assert annotated.loc[0, "pearson_pvalue"] < 0.05
