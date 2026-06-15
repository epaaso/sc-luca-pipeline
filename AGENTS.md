# Agent Notes For This Repo

This repository contains the unified Nextflow pipeline for the LUCA single-cell
analysis. Read this file before changing code or launching jobs.

## Current Objective And Status

The implemented workflow phases are:

- `RAYTUNE`: scVI hyperparameter search.
- `ATLAS`: reference scVI/scANVI training.
- `SURGERY`: annotation of one query dataset.
- `SUBCLUSTER`: early/late tumor-epithelial clustering.
- `GRAPH`: early/late global ARACNe, sample ecotype clustering, and
  ecotype-specific ARACNe with Pearson annotations.
- `PIPELINE`: chains the model, cohort surgery, subclustering, and graph phases.

The graph implementation is based on
`../sc-luca-explore/nb_graphGeneration` and
`../sc-luca-explore/nb_ecotypeClustering`.

As of June 14, 2026:

- Unit and CPU integration suite passes: `16 passed`.
- Standalone `GRAPH` was validated through Slurm with a synthetic cohort.
- Standalone graph checkpoint execution was validated.
- Remote ARACNe and Java runtimes were provisioned.
- A real biological graph run has **not** been completed.
- The complete `PIPELINE` has not been validated end to end after adding
  `GRAPH`. Its last stub attempt exposed an existing remote staging problem in
  `PREPARE_DATASET`: the staged download script could not be found in the
  remote work directory.

## Ground Rules

- Do not reset, delete, or overwrite user outputs unless explicitly requested.
- The worktree is intentionally dirty. Preserve all existing changes and
  untracked implementation files.
- The graph implementation is not committed. Relevant changed/untracked files
  include `main.nf`, `nextflow.config`, `bin/run_subcluster.py`,
  `bin/run_graph_phase.py`, `configs/graph_default.yaml`, and graph/subcluster
  tests and documentation.
- Prefer small, reversible edits.
- Use `rg`, `find`, and `sed` for inspection and `apply_patch` for manual edits.
- Before restarting long jobs, check for reusable models and use `-resume`.
- Do not describe a synthetic smoke run as a biological result.
- ARACNe is CPU/Java work. Under `remote_gpu`, it runs through Slurm but does
  not request a GPU.

## Repository And Reference Paths

```text
Pipeline root:
/datos/home/epaaso/REPOS/sc-luca-pipeline

Reference notebooks and historical outputs:
/datos/home/epaaso/REPOS/sc-luca-explore
/datos/migccl/neto_maestria/luca_explore

Nextflow Slurm work directory:
/datos/home/epaaso/slurm-gpu-jobs/nf-work
```

Remote runtime paths:

```text
Training container:
/data/containers/scvi-raytune-py313-cu12.sif

Subcluster/graph container:
/data/containers/sc-luca-subcluster-cu12.sif

Shared Python dependency layer:
/data/containers/sc-luca-subcluster-py311
/data/containers/sc-luca-subcluster-libs

ARACNe:
/data/containers/ARACNe-AP/dist/aracne.jar

Java 11:
/data/containers/java/bin/java
```

Remote host used for provisioning and checks:

```bash
ssh epaaso@inmegen-gpu
```

## Important Implementation Files

- `main.nf`: processes, named workflows, and default chained `PIPELINE`.
- `nextflow.config`: `local`, `local_cpu`, `local_gpu`, and `remote_gpu`.
- `bin/run_subcluster.py`: cohort normalization, clustering, diagnostics, and
  mapping application.
- `bin/run_graph_phase.py`: graph preparation, ARACNe execution, ecotype
  clustering, and Pearson annotation.
- `configs/subcluster_default.yaml`: standalone subcluster defaults.
- `configs/subcluster_cohort_default.json`: canonical atlas/query manifest.
- `configs/graph_default.yaml`: standalone graph defaults.
- `tests/test_subcluster*.py` and `tests/test_graph_phase.py`: current tests.

## Data Contracts Between Phases

`SUBCLUSTER` raw output includes:

```text
subcluster/query_latent.h5ad
subcluster/query_latent_tumor_early.h5ad
subcluster/query_latent_tumor_late.h5ad
```

The mapping process output always contains:

```text
mapping_results/query_latent_adjusted.h5ad
```

It is also copied into the experiment shared output as:

```text
<shared_run_root>/<experiment_name>/subcluster/query_latent_adjusted.h5ad
```

Its `obs` contains `cell_type_adjusted`. When no mapping is supplied, adjusted
labels equal the base `cell_type` labels. This H5AD is the canonical chained
input to `GRAPH`.

Standalone `GRAPH` accepts an H5AD or CSV with configurable columns. Defaults:

```text
cell type: cell_type_adjusted, with fallback to cell_type
sample:    batch
dataset:   dataset
stage:     stage
```

Stages are split strictly into early `I`/`II` and late
`III`/`III or IV`/`IV`. Unknown stages are excluded and counted in QC.

## Execution Commands

Run all tests:

```bash
pytest -q
```

Standalone remote FAISS subclustering:

```bash
nextflow run main.nf -entry SUBCLUSTER \
  -params-file configs/subcluster_default.yaml \
  -profile remote_gpu -resume
```

Standalone complete graph phase:

```bash
nextflow run main.nf -entry GRAPH \
  -params-file configs/graph_default.yaml \
  -profile remote_gpu \
  --graph_input /path/to/query_latent_adjusted.h5ad \
  --experiment_name meaningful_graph_run \
  -resume
```

Always override `--graph_input` for meaningful runs. The path currently stored
in `configs/graph_default.yaml` points to a small historical subcluster smoke
artifact and exists only as an executable example.

Graph checkpoints use `--graph_step`:

```text
prepare
global-aracne
ecotype
ecotype-aracne
all
```

For steps after `prepare`, pass `--graph_global_dir`. For
`ecotype-aracne`, also pass `--ecotype_dir`.

Default full graph settings are production-oriented:

```text
aracne_bootstraps: 500
aracne_pvalue: 1E-8
ecotype_k: 4
ecotype_resolution: 0.2
```

Use `--aracne_dry_run true` to validate matrices and ecotypes without Java.

## Validated Artifacts

Synthetic input used for remote validation:

```text
/datos/home/epaaso/slurm-gpu-jobs/graph_inputs/synthetic_graph.h5ad
```

Successful final remote graph smoke:

```text
/datos/home/epaaso/slurm-gpu-jobs/graph/runs/graph_remote_final_20260613/graph
```

The smoke used `3` ARACNe bootstraps and `pvalue=1E-2`. It recovered two
six-sample ecotypes in each stage and produced all expected files. Its final
networks contain zero significant edges; that is expected for this small
synthetic validation and is not a production result.

Latest available subcluster smoke:

```text
/datos/home/epaaso/slurm-gpu-jobs/subcluster/runs/faiss_remote_gpu_run_20260609_v3/subcluster
```

That smoke has only one early sample and one late sample, so it is not suitable
for a meaningful ecotype/ARACNe run.

## Known Issues And Next Work

1. Run `SUBCLUSTER` on the actual complete atlas plus Bishoff, Deng, Hu, and
   Zuani cohort, curate complete early/late mappings, and generate
   `query_latent_adjusted.h5ad`.
2. Run a production `GRAPH` job with the default `500` bootstraps.
3. Diagnose and fix the complete `PIPELINE` remote staging failure in
   `PREPARE_DATASET`, then validate the full chain with `-stub-run` and a small
   real configuration.
4. Current `GRAPH` execution runs all early/late and ecotype ARACNe work inside
   one Nextflow process. A future scalability improvement is to fan out ARACNe
   jobs by stage and ecotype.
5. ARACNe bootstrap files are retained and can consume substantial storage in
   production runs. Add an explicit retention/cleanup option before changing
   this behavior.
6. The next biological downstream phase should run differential expression and
   CellphoneDB separately per patient ecotype, rather than calculating global
   markers and filtering afterward.

## Final Outputs

- `ATLAS`: `ref_latent.h5ad`, model directories, metrics, histories.
- `SURGERY`: query H5ADs and prediction CSVs.
- `RAYTUNE`: best configuration/result and trial summary.
- `SUBCLUSTER`: merged/stage H5ADs, Leiden/UMAP CSVs, diagnostics, mappings.
- `GRAPH`: global matrices/networks, sample ecotype memberships, abundance
  matrices, sample kNN graphs, and ecotype-specific MI/Pearson networks.
