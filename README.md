# LUCA Single-Cell Coabundance Networks Pipeline

The main objective of this project is to perform an ecological analysis of the cell types in NSCLC lung tumor tissues from scRNA-seq data, in both early and late tumor stages. We obtain coabundance metrics, specifically relying on Mutual Information (MI), as it captures non-linear correlations in the data, and the methodology in the [ARACNE-AP](https://github.com/califano-lab/ARACNe-AP) package ensures a smaller possibility of spurious correlations.

This repository serves as the unified Nextflow automation of preprocessing, model training, annotation, tumor/epithelial subclustering, and downstream analysis. It replaces the old Jupyter notebook workflow with a robust, scalable pipeline.

![scRNA-seq Workflow](Workflow_dark.png "scRNA-seq Workflow")


## Pipeline Execution Modes

The pipeline supports two execution models: the default auto-chained end-to-end workflow, or individual standalone workflows.

### 1. End-to-End Auto-Chained Pipeline (`PIPELINE` Mode)

This is the default execution workflow. It chains all stages together:
1. **Hyperparameter Search**: Runs hyperparameter optimization via Ray Tune over the SCVI model parameters.
2. **Config Extraction**: Identifies the best trial and exports the optimal settings.
3. **Atlas Training**: Merges the optimized parameters into the baseline configuration and trains a reference SCANVI model.
4. **Cohort Dataset Surgery**: Maps the configured Bishoff, Deng, Hu, and Zuani datasets onto the newly trained SCANVI atlas.
5. **Subclustering**: Merges atlas/query latent representations, splits early and late disease, and runs tumor/epithelial FAISS-Leiden subclustering.
6. **Graph And Ecotype Analysis**: Builds early/late global ARACNe networks, clusters samples into ecotypes, then builds and Pearson-annotates one ARACNe network per ecotype.

To run the complete chained pipeline:
```bash
# Run end-to-end pipeline on Slurm GPU cluster
nextflow run main.nf -profile remote_gpu
```

The chained pipeline is primarily intended for `remote_gpu`, because model training and the default subclustering backend require a GPU.

> **Validation status:** standalone subclustering and graph execution have been
> validated, including a real Slurm graph smoke run on synthetic data. A
> production biological graph run and a complete end-to-end `PIPELINE` run are
> still pending. See [Current Development Status](#current-development-status).

### 2. Standalone Workflows

You can also run individual phases of the pipeline separately by using the `-entry` flag:

* **Hyperparameter Tuning (`RAYTUNE` Workflow)**: Explores SCVI model parameter search spaces via Ray Tune.
  ```bash
  nextflow run main.nf -entry RAYTUNE -params-file configs/raytune_default.yaml -profile remote_gpu
  ```

* **Atlas Generation (`ATLAS` Workflow)**: Trains the scVI/scANVI reference models to embed and cluster reference datasets.
  ```bash
  nextflow run main.nf -entry ATLAS -params-file configs/atlas_default.yaml -profile remote_gpu
  ```

* **Dataset Annotation (`SURGERY` Workflow)**: Annotates external query datasets via architectural surgery using a pre-trained reference model.
  ```bash
  nextflow run main.nf -entry SURGERY -params-file configs/surgery_default.yaml -profile local
  ```

* **Tumor/Epithelial Subclustering (`SUBCLUSTER` Workflow)**: Runs independently from a cohort manifest or a canonical premerged latent H5AD.
  ```bash
  # Local CPU Scanpy
  nextflow run main.nf -entry SUBCLUSTER \
    -params-file configs/subcluster_default.yaml \
    -profile local_cpu

  # Local GPU FAISS
  nextflow run main.nf -entry SUBCLUSTER \
    -params-file configs/subcluster_default.yaml \
    -profile local_gpu

  # Slurm GPU FAISS
  nextflow run main.nf -entry SUBCLUSTER \
    -params-file configs/subcluster_default.yaml \
    -profile remote_gpu
  ```

* **Graph And Ecotype Analysis (`GRAPH` Workflow)**: Runs the complete graph phase from a merged H5AD or cell-level CSV.
  ```bash
  nextflow run main.nf -entry GRAPH \
    -params-file configs/graph_default.yaml \
    -profile remote_gpu \
    --graph_input /path/to/subcluster/query_latent_adjusted.h5ad
  ```

  Run individual checkpoints with `--graph_step prepare`, `global-aracne`, `ecotype`, or `ecotype-aracne`. For steps after preparation, set `graph_global_dir`; for ecotype-specific ARACNe, also set `ecotype_dir` to prior outputs. `--graph_step all` is the default.

## Dataset Catalog (`metadata/dsets.csv`)

The datasets used in this pipeline are standardized and documented under [metadata/dsets.csv](metadata/dsets.csv). This CSV acts as the primary registry mapping quality control thresholds, platform chemistries, biological/clinical details, and publication DOIs for all integrated single-cell RNA-seq studies.

Key columns in the registry include:
- `id`: The unique name identifier of the dataset (e.g., `Chen_Zhang_2020_NSCLC`, `Kim_Lee_2020_LUAD`, `Zuani_2024_NSCLC`).
- `input_adata`: The expected local path to the raw/processed AnnData file.
- `min_counts`, `max_counts`, `min_genes`, `max_genes`, `max_pct_mito`: Quality control (QC) thresholds applied during preprocessing.
- `integrated cells n`: Total cells included from the study.
- `doi`, `journal`, `names`: Citation and publication registry information.
- `stage`, `disease`, `cell_sorting`, `prior_treatment`: Biological details and clinical variables of the dataset.

---

### Key Execution Parameters

When running the pipeline or individual workflows, you can override default values in the configuration files or pass them dynamically on the command line (e.g., `--input_h5ad "/path/to/file"`). The most critical parameters include:

#### 1. Input Datasets & Paths
- `--input_h5ad`: Path to the input atlas/training AnnData file (e.g., `/data/luca_atlas/extended.h5ad` or `extended_tumor_hvg.h5ad`). If the file does not exist, the pipeline automatically retrieves and prepares it.
- `--query_h5ad` (in `PIPELINE` mode) or `--input_h5ad` (in `SURGERY` mode): Path to the query AnnData file to project.
- `--reference_model` (in `SURGERY` standalone mode): Folder containing the pre-trained SCANVI reference model (`model.pt`).
- `--shared_run_root`: Root path where all run logs, TensorBoard files, checkpoints, and final `.h5ad` files are stored.
- `--experiment_name`: Unique name identifier for the current run, used to structure output directories.

#### 2. SCVI/SCANVI Model settings
- `gene_likelihood` (under `scvi_model_params`): Set to `nb` (Negative Binomial) to model biological zeros. Avoid using `zinb` for this atlas context.
- `lr` (under `scvi_plan_kwargs` and `scanvi_plan_kwargs`): Set to `1e-4` to prevent numerical training instabilities (NaNs).
- `gradient_clip_val` (under `scvi_trainer_kwargs` and `scanvi_trainer_kwargs`): Set to `10.0` to preserve training stability.

#### 3. Scaling & Multiprocessing Settings
- `--dl_num_workers`: Set to `2` to prevent PyTorch dataloaders from triggering "Too many open files" errors.
- `torch_multiprocessing_sharing_strategy`: Set to `file_system` to resolve resource sharing limits during dataloader multi-worker training.

---

## Subclustering

The `SUBCLUSTER` workflow automates the core behavior from `Tumor_subcluster_wu.ipynb`:

1. Load the atlas latent coordinates from `ref_latent.h5ad:X`.
2. Load surgery coordinates from `query_latent.h5ad:obsm["X_scVI"]`.
3. Normalize cell type, dataset, batch, and stage metadata.
4. Select predicted tumor cells plus the configured epithelial populations.
5. Split cells strictly into early (`I`, `II`) and late (`III`, `IV`, `III or IV`) cohorts.
6. Run independent FAISS or Scanpy neighbor, Leiden, and UMAP calculations.

Missing or unrecognized stages are excluded and reported in `qc_summary.json`.

### Input Contracts

Set exactly one of:

- `subcluster_cohort_manifest`: JSON manifest describing an atlas and query latent artifacts. The canonical example is [configs/subcluster_cohort_default.json](configs/subcluster_cohort_default.json).
- `subcluster_input_h5ad`: Premerged compact latent H5AD whose `X` is the latent matrix and whose `obs` contains `cell_type`, `dataset`, `batch`, and `stage`.

Manifest query entries may provide `obs_csv` and `obs_join: position` or `index` for legacy latent H5AD files without observation metadata. Stage fields are configured per dataset; the default cohort uses `Pathological stage` for Deng, `Clinical Stage` for Hu, `stage` for Zuani, and an explicit `III or IV` override for Bishoff.

### Engines And Profiles

- `local_cpu` selects `subcluster_engine=scanpy` and `subcluster_backend=cpu`, then runs the subclustering container with Docker.
- `local_gpu` selects `subcluster_engine=faiss` and `subcluster_backend=gpu`, then runs Docker with `--gpus all`. The host must provide NVIDIA Container Toolkit.
- `remote_gpu` selects `subcluster_engine=faiss` and `subcluster_backend=gpu`, submits to Slurm, and runs the configured GPU SIF with Singularity `--nv`.
- GPU mode runs FAISS neighbor search on GPU and defaults to `subcluster_gpu_postprocess=cpu` for Leiden and UMAP. Set it to `rapids` only when the container CUDA toolkit supports the target GPU architecture.
- FAISS uses scalable IVF parameters derived from cohort size unless `subcluster_faiss_nlist` or `subcluster_faiss_nprobe` is set.
- Scanpy is the CPU reference implementation and rejects GPU execution. Engine and backend parameters remain command-line overridable for advanced runs.

Shared scientific defaults are 30 neighbors, Leiden resolution `0.5`, UMAP `min_dist=0.03`, UMAP spread `0.6`, and seed `42`.

### Mapping Checkpoint

Raw cluster IDs are always emitted. Semantic names are applied only when complete per-stage mapping JSON files are supplied through `early_mapping_json` and `late_mapping_json`:

```json
{
  "0": "Tumor AT2-like",
  "1": "Tumor EMT"
}
```

Mappings must contain exactly the observed cluster IDs. Without mappings, the run completes successfully and skips adjusted per-dataset prediction exports. Re-run with `-resume` after curation to reuse cached clustering.

### Subclustering Outputs

Each experiment writes to `<shared_run_root>/<experiment_name>/subcluster/`:

- `query_latent.h5ad`: compact normalized merged cohort.
- `query_latent_tumor_early.h5ad` and `query_latent_tumor_late.h5ad`: clustered stage subsets.
- `atlas_<stage>_leiden.csv`, `atlas_<stage>_umap.csv`, and `atlas_<stage>_uparams.json`.
- `qc_summary.json`, `versions.json`, contingency CSVs, UMAP plots, and heatmaps.
- `<dataset>_predicted_leiden_<stage>.csv` when a valid mapping is supplied.
- `query_latent_adjusted.h5ad`: merged cohort with `cell_type_adjusted`; base labels are retained for stages without mappings. Inside the Nextflow mapping task this file is under `mapping_results/`, then it is copied into the experiment's shared `subcluster/` output directory.

---

## Graph Generation And Ecotypes

The `GRAPH` workflow replaces the early/late graph-generation and abundance-clustering notebooks. It accepts a cell-level CSV or H5AD and:

1. Counts `cell_type_adjusted` per sample, keeping type/sample combinations with at least two cells.
2. Writes early (`I`, `II`) and late (`III`, `III or IV`, `IV`) sample-by-cell-type matrices.
3. Runs global ARACNe-AP independently for early and late disease.
4. Normalizes sample abundances, builds a sample kNN graph, and discovers ecotypes with Leiden.
5. Subsets each stage matrix by ecotype, runs ARACNe-AP again, and annotates resulting MI edges with Pearson correlation and p-value.

The default graph input columns match `query_latent_adjusted.h5ad`: `cell_type_adjusted`, `batch`, `dataset`, and `stage`. Override `graph_cell_type_key`, `graph_sample_key`, `graph_dataset_key`, or `graph_stage_key` for another schema. When `cell_type_adjusted` is absent, the driver falls back to `cell_type`.

ARACNe is CPU/Java software. The `remote_gpu` profile still submits graph work through Slurm, but graph tasks do not request a GPU. The remote runtime expects:

- `/data/containers/ARACNe-AP/dist/aracne.jar`
- `/data/containers/java/bin/java`
- `/data/containers/sc-luca-subcluster-cu12.sif` and its shared Python dependency layer

Important controls are `aracne_bootstraps` (default `500`), `aracne_pvalue` (default `1E-8`), `ecotype_k` (default `4`), and `ecotype_resolution` (default `0.2`). Use `aracne_dry_run: true` to validate matrix and ecotype preparation without launching Java.

Outputs are written to `<shared_run_root>/<experiment_name>/graph/`:

- `global/groups_{early,late}.csv`, count matrices, QC, and global ARACNe networks.
- `ecotype/membership_{early,late}.csv`, normalized abundance matrices, sample kNN GraphML files, and cluster plots.
- `ecotype_graphs/{early,late}/cluster_XX/` with ARACNe inputs, networks, QC, and `*_MI_pearson.txt` annotations.

### Running From Subcluster Results

After cluster mappings are reviewed, use the merged adjusted H5AD as the graph
input:

```bash
nextflow run main.nf -entry GRAPH \
  -params-file configs/graph_default.yaml \
  -profile remote_gpu \
  --graph_input /path/to/subcluster/query_latent_adjusted.h5ad \
  --experiment_name luca_graph_production \
  -resume
```

The default `500` ARACNe bootstraps are intended for production. For a quick
runtime smoke test, override `--aracne_bootstraps 3`; do not interpret such a
smoke run biologically.

Always override `graph_input` for a meaningful run. The input currently present
in `configs/graph_default.yaml` is a small historical smoke artifact.

### Graph Checkpoint Examples

Prepare matrices without running ARACNe:

```bash
nextflow run main.nf -entry GRAPH \
  -params-file configs/graph_default.yaml \
  -profile remote_gpu \
  --graph_input /path/to/query_latent_adjusted.h5ad \
  --graph_step prepare \
  --aracne_dry_run true
```

Re-run ecotype discovery from prepared global matrices:

```bash
nextflow run main.nf -entry GRAPH \
  -params-file configs/graph_default.yaml \
  -profile remote_gpu \
  --graph_input /path/to/query_latent_adjusted.h5ad \
  --graph_step ecotype \
  --graph_global_dir /path/to/prior/graph/global
```

Run per-ecotype ARACNe from reviewed memberships:

```bash
nextflow run main.nf -entry GRAPH \
  -params-file configs/graph_default.yaml \
  -profile remote_gpu \
  --graph_input /path/to/query_latent_adjusted.h5ad \
  --graph_step ecotype-aracne \
  --graph_global_dir /path/to/prior/graph/global \
  --ecotype_dir /path/to/prior/graph/ecotype
```

### Container Build

Build the reproducible CPU/GPU image:

```bash
docker build -t sc-luca-subcluster:latest -f containers/subcluster/Dockerfile .
singularity build /data/containers/sc-luca-subcluster-cu12.sif containers/subcluster/subcluster.def
```

For Slurm/Singularity execution, publish or convert that image to the path configured by `subcluster_container`, defaulting to `/data/containers/sc-luca-subcluster-cu12.sif` under `remote_gpu`.

On the current GPU host, `remote_gpu` uses the FAISS/RAPIDS base SIF plus a small shared Python dependency layer:

```bash
singularity pull /data/containers/faiss-cugraph-24-12.sif docker://netopaas/faiss:cugraph-24-12
mkdir -p /data/containers/sc-luca-subcluster-py311
singularity exec --bind /data /data/containers/faiss-cugraph-24-12.sif \
  /opt/conda/bin/python -m pip install --target /data/containers/sc-luca-subcluster-py311 \
  igraph==0.11.9 leidenalg==0.10.2
ln -s /data/containers/faiss-cugraph-24-12.sif /data/containers/sc-luca-subcluster-cu12.sif
```

The current local host does not expose an NVIDIA Docker runtime, so `local_gpu` must be validated on a host configured with NVIDIA Container Toolkit.

The remote graph runtime also requires:

```text
/data/containers/ARACNe-AP/dist/aracne.jar
/data/containers/java/bin/java
```

---

## Automatic Dataset Downloading & Preprocessing

To make execution seamless and remote-friendly (especially on environments without access to local shared mounts like `/datos`), the pipeline automatically downloads and prepares input files if they do not exist locally:

1. **LUCA Reference Atlas**:
   - If the files `/data/luca_atlas/extended.h5ad` or `/data/luca_atlas/extended_tumor_hvg.h5ad` are missing, Nextflow will trigger `PREPARE_DATASET` for the `extended` or `extended_tumor_hvg` dataset.
   - The pipeline retrieves the full extended atlas directly from CELLxGENE's stable dataset asset portal (`https://datasets.cellxgene.cziscience.com/173984ce-d33a-46b6-ae96-4be47f6c67e8.h5ad`).
   - For `extended_tumor_hvg`, the script automatically subsets the downloaded atlas in memory to `tumor_primary` cells, filters for specific UICC stages (`I` through `IV`), excludes conflicting studies, filters for highly variable genes, and prepares raw counts.

2. **Query/Surgery Datasets**:
   - If the input query file specified for dataset surgery (e.g., `Bishoff`, `Hu`, `Zuani`, `Deng`) does not exist at its designated location, the pipeline automatically fetches the raw matrices/barcodes/features from their public endpoints (EBI BioStudies, Figshare, GEO, CodeOcean) and runs the corresponding filtering pipeline defined in `bin/download_and_preprocess_dataset.py` before executing the surgery mapping.

---

## Folder Structure

- `main.nf` - The primary Nextflow workflow definition containing `ATLAS`, `SURGERY`, `RAYTUNE`, `SUBCLUSTER`, `GRAPH`, and `PIPELINE`.
- `nextflow.config` - Defines `local`, `local_cpu`, `local_gpu`, and `remote_gpu` executor/container profiles.
- `bin/run_subcluster.py` - Cohort preparation, FAISS/Scanpy subclustering, diagnostics, and mappings.
- `bin/run_graph_phase.py` - Global/ecotype matrix preparation, ARACNe, Leiden ecotypes, and Pearson annotations.
- `configs/` - Base YAML parameter files, graph defaults, and the canonical subclustering cohort manifest.
- `containers/subcluster/` - Reproducible FAISS/RAPIDS subclustering image recipe.
- `tests/` - Unit and CPU integration tests for subclustering and graph preparation.
- `metadata/` - Contains dataset catalog files like `dsets.csv`.
- `AGENTS.md` and `SKILLS.md` - Context for interacting with the AI agent system to develop or monitor the pipeline.

---

## Current Development Status

Validated as of June 14, 2026:

- `pytest -q`: `16 passed`.
- Remote Slurm FAISS subclustering smoke completed previously.
- Standalone remote `GRAPH` and graph checkpoint execution completed.
- The final synthetic graph smoke is under:
  `/datos/home/epaaso/slurm-gpu-jobs/graph/runs/graph_remote_final_20260613/graph`.

The synthetic graph smoke used two six-sample ecotypes per stage, three ARACNe
bootstraps, and a relaxed p-value. All expected artifacts were emitted, but the
networks contain no significant edges. This validates execution only.

Still pending:

- Run subclustering on the complete biological cohort and curate early/late
  cluster mappings.
- Run production graph generation with the adjusted biological labels and
  default `500` bootstraps.
- Validate the complete chained `PIPELINE`. The latest stub attempt failed
  before graph execution because `PREPARE_DATASET` could not find its staged
  download script in the remote Slurm work directory.
- Scale ARACNe by fanning stage/ecotype jobs into separate Nextflow tasks if
  production runtime requires it.
