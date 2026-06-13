import importlib.util
import json
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest


MODULE_PATH = pathlib.Path(__file__).parents[1] / "bin" / "run_subcluster.py"
SPEC = importlib.util.spec_from_file_location("run_subcluster_integration", MODULE_PATH)
subcluster = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(subcluster)


def make_premerged(path: pathlib.Path) -> None:
    rng = np.random.default_rng(42)
    early = np.vstack([rng.normal(-2, 0.3, (30, 5)), rng.normal(2, 0.3, (30, 5))])
    late = np.vstack([rng.normal(-1, 0.3, (30, 5)), rng.normal(3, 0.3, (30, 5))])
    matrix = np.vstack([early, late]).astype(np.float32)
    obs = pd.DataFrame(
        {
            "cell_type": ["Tumor test"] * 120,
            "dataset": ["D1"] * 60 + ["D2"] * 60,
            "batch": ["b1"] * 60 + ["b2"] * 60,
            "stage": ["I"] * 60 + ["III"] * 60,
        },
        index=[f"cell_{i}" for i in range(120)],
    )
    ad.AnnData(X=matrix, obs=obs).write_h5ad(path)


@pytest.mark.parametrize("engine", ["scanpy", "faiss"])
def test_cpu_subcluster_end_to_end(tmp_path, engine):
    required = ["igraph", "leidenalg"]
    if engine == "faiss":
        required.extend(["faiss", "umap"])
    for module in required:
        if importlib.util.find_spec(module) is None:
            pytest.skip(f"{module} is not installed")

    input_path = tmp_path / "cohort.h5ad"
    make_premerged(input_path)
    output_dir = tmp_path / engine
    cfg = {
        "subcluster_input_h5ad": str(input_path),
        "subcluster_cohort_manifest": None,
        "subcluster_output_dir": str(output_dir),
        "subcluster_engine": engine,
        "subcluster_backend": "cpu",
        "subcluster_seed": 42,
        "subcluster_n_neighbors": 8,
        "subcluster_leiden_resolution": 0.5,
        "subcluster_umap_min_dist": 0.03,
        "subcluster_umap_spread": 0.6,
    }
    subcluster.run(cfg)

    assert (output_dir / "query_latent.h5ad").exists()
    assert (output_dir / "query_latent_tumor_early.h5ad").exists()
    assert (output_dir / "atlas_late_leiden.csv").exists()
    with (output_dir / "qc_summary.json").open() as handle:
        qc = json.load(handle)
    assert qc["stages"]["early"]["n_obs"] == 60
    assert qc["stages"]["late"]["n_obs"] == 60
