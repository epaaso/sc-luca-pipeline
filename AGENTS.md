# Agent Notes For This Repo

This repo contains the unified Nextflow pipeline for the LUCA single-cell analysis.
The workflows handle the deep learning portions of the project: `ATLAS` (model training), `SURGERY` (dataset annotation), and `RAYTUNE` (hyperparameter search).

## Ground Rules

- Do not reset, delete, or overwrite user outputs unless the user explicitly asks.
- Prefer small, reversible edits in the scripts and configurations.
- Use `rg`/`find`/`sed` for inspection and `apply_patch` for manual file edits.
- Before restarting long GPU jobs, check whether a reusable model already exists.

## Important Paths

Unified pipeline root:

```text
/home/epaaso/REPOS/sc-luca-pipeline
```

Shared LUCA working tree:

```text
/datos/migccl/neto_maestria/luca_explore
```

Container image commonly used:

```text
/data/containers/scvi-raytune-py313-cu12.sif
```

## Pipeline Execution

The pipeline is implemented in Nextflow, using `main.nf` with multiple named entry workflows (`ATLAS`, `SURGERY`, `RAYTUNE`).
Executions can be parameterized with external YAML configurations.
The pipeline supports `local` and `remote_gpu` profiles via `nextflow.config` to gracefully handle different compute environments.

## Final Outputs

Outputs depend on the workflow run:
- **ATLAS**: `ref_latent.h5ad`, `scvi_model/`, `scanvi_model/`, and `metrics.json`.
- **SURGERY**: `query.h5ad`, `query_latent.h5ad`, `*_predicted.csv`, and mapped `surgery_model/`.
- **RAYTUNE**: `best_config.json`, `best_result.json`, `summary.csv`.

## Architectural Insights & Future Work

### Per-Cluster CellphoneDB & Differential Expression
**Current State**:
- `modal_DE.py` computes differential expression globally across all patients in a stage.
- `run_cellphone_full.py` runs CellphoneDB globally and then filters results based on cell types present in an ARACNE patient cluster.

**The Flaw**:
ARACNE clusters are subsets of patients. Using global marker genes assumes cell types (e.g., Macrophages) have identical expression profiles across all patient groups, ignoring intra-cluster biology.

**Next Steps / Required Fix**:
To properly capture intra-cluster interactions:
1. Divide single-cell data by patient ARACNE clusters.
2. Run differential expression (`modal_DE.py` / `aggregate_markers.py`) *separately* for each cluster.
3. Run CellphoneDB separately per cluster using these cluster-specific marker gene profiles.
