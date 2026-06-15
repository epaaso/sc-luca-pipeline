import importlib.util
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest


MODULE_PATH = pathlib.Path(__file__).parents[1] / "bin" / "run_subcluster.py"
SPEC = importlib.util.spec_from_file_location("run_subcluster", MODULE_PATH)
subcluster = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(subcluster)


def test_canonicalize_stage():
    assert subcluster.canonicalize_stage("IA1") == "I"
    assert subcluster.canonicalize_stage("IIB") == "II"
    assert subcluster.canonicalize_stage("IIIA") == "III"
    assert subcluster.canonicalize_stage("III or IV") == "III or IV"
    assert subcluster.canonicalize_stage("unknown") is None


def test_cell_selection_and_stage_split():
    obs = pd.DataFrame(
        {
            "cell_type": ["Tumor LUAD", "Club", "T cell CD4", "Ciliated"],
            "stage": ["I", "III", None, "bad"],
        }
    )
    assert subcluster.select_subcluster_cells(obs).tolist() == [True, True, False, True]
    masks = subcluster.split_stage_masks(obs)
    assert masks["early"].tolist() == [True, False, False, False]
    assert masks["late"].tolist() == [False, True, False, False]
    assert masks["unknown"].tolist() == [False, False, True, True]


def test_effective_faiss_params_are_bounded():
    assert subcluster.effective_faiss_params(100, 200) == (99, 10, 10)
    assert subcluster.effective_faiss_params(100, 30, nlist=200, nprobe=300) == (30, 100, 100)


def test_symmetric_affinity_graph_is_symmetric_and_positive():
    indices = np.array([[1, 2], [0, 2], [0, 1]])
    distances = np.array([[1.0, 2.0], [1.0, 4.0], [2.0, 4.0]])
    graph = subcluster.symmetric_affinity_graph(indices, distances, 3)
    assert graph.shape == (3, 3)
    assert np.allclose(graph.toarray(), graph.toarray().T)
    assert np.all(graph.data > 0)


def test_mapping_validation_is_strict():
    assert subcluster.validate_mapping(["0", "1"], {"0": "A", "1": "B"}) == {"0": "A", "1": "B"}
    with pytest.raises(ValueError, match="missing clusters"):
        subcluster.validate_mapping(["0", "1"], {"0": "A"})


def test_apply_mapping_preserves_nonselected_and_source_indexes(tmp_path):
    obs = pd.DataFrame(
        {
            "source_index": ["a", "b", "c"],
            "cell_type": ["Tumor old", "T cell", "Tumor old"],
            "dataset": ["D1", "D1", "D2"],
            "batch": ["b1", "b1", "b2"],
            "stage": ["I", "I", "III"],
        },
        index=["D1::a", "D1::b", "D2::c"],
    )
    merged = ad.AnnData(X=np.ones((3, 2)), obs=obs)
    subcluster.apply_mapping_exports(merged, "early", ["D1::a"], ["0"], {"0": "Tumor new"}, tmp_path)
    d1 = pd.read_csv(tmp_path / "D1_predicted_leiden_early.csv")
    assert d1["index"].tolist() == ["a", "b"]
    assert d1["cell_type_adjusted"].tolist() == ["Tumor new", "T cell"]


def test_premerged_validation_normalizes_stages():
    obs = pd.DataFrame(
        {
            "cell_type": ["Tumor A", "Tumor B"],
            "dataset": ["D", "D"],
            "batch": ["b", "b"],
            "stage": ["IA1", "IIIA"],
        },
        index=["a", "b"],
    )
    result = subcluster.validate_premerged(ad.AnnData(X=np.ones((2, 3)), obs=obs))
    assert result.obs["stage"].tolist() == ["I", "III"]
    assert result.obs["source_index"].tolist() == ["a", "b"]


def test_compatible_h5ad_writer_normalizes_string_arrays(tmp_path):
    obs = pd.DataFrame(
        {"cell_type": pd.array(["Tumor A", "Tumor B"], dtype="string")},
        index=pd.Index(pd.array(["a", "b"], dtype="string")),
    )
    output = tmp_path / "compatible.h5ad"
    subcluster.write_compatible_h5ad(ad.AnnData(X=np.ones((2, 3)), obs=obs), output)
    assert ad.read_h5ad(output).obs_names.tolist() == ["a", "b"]


def test_manifest_entry_supports_external_obs_position_join(tmp_path):
    latent_path = tmp_path / "latent.h5ad"
    obs_path = tmp_path / "obs.csv"
    ad.AnnData(X=np.ones((2, 4)), obs=pd.DataFrame(index=["0", "1"])).write_h5ad(latent_path)
    pd.DataFrame(
        {
            "index": ["D_a", "D_b"],
            "cell_type_tumor": ["Tumor A", "Club"],
            "sample": ["s1", "s2"],
            "stage": ["I", "III"],
        }
    ).to_csv(obs_path, index=False)
    result = subcluster.load_manifest_entry(
        {
            "name": "D",
            "latent_h5ad": str(latent_path),
            "obs_csv": str(obs_path),
            "obs_csv_index_col": None,
            "obs_join": "position",
            "source_index_key": "index",
            "cell_type_key": "cell_type_tumor",
            "batch_key": "sample",
            "stage_key": "stage",
        }
    )
    assert result.obs["source_index"].tolist() == ["D_a", "D_b"]
    assert result.obs["stage"].tolist() == ["I", "III"]


def test_apply_mappings_from_raw_results(tmp_path):
    raw_dir = tmp_path / "raw"
    output_dir = tmp_path / "mapped"
    raw_dir.mkdir()
    obs = pd.DataFrame(
        {
            "source_index": ["a", "b"],
            "cell_type": ["Tumor old", "T cell"],
            "dataset": ["D", "D"],
            "batch": ["b", "b"],
            "stage": ["I", "I"],
        },
        index=["D::a", "D::b"],
    )
    merged = ad.AnnData(X=np.ones((2, 2)), obs=obs)
    merged.write_h5ad(raw_dir / "query_latent.h5ad")
    stage = merged[["D::a"]].copy()
    stage.obs["leiden"] = pd.Categorical(["0"])
    stage.write_h5ad(raw_dir / "query_latent_tumor_early.h5ad")
    mapping_path = tmp_path / "early.json"
    mapping_path.write_text('{"0": "Tumor new"}')

    subcluster.apply_mappings(
        raw_dir,
        {"early_mapping_json": str(mapping_path), "late_mapping_json": None},
        output_dir,
    )
    mapped = pd.read_csv(output_dir / "D_predicted_leiden_early.csv")
    assert mapped["cell_type_adjusted"].tolist() == ["Tumor new", "T cell"]
    adjusted = ad.read_h5ad(output_dir / "query_latent_adjusted.h5ad")
    assert adjusted.obs["cell_type_adjusted"].tolist() == ["Tumor new", "T cell"]
