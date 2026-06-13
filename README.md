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

To run the complete chained pipeline:
```bash
# Run end-to-end pipeline on Slurm GPU cluster
nextflow run main.nf -profile remote_gpu
```

The chained pipeline is primarily intended for `remote_gpu`, because model training and the default subclustering backend require a GPU.

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

*Note: The remaining phases—mutual information extraction, network generation with ARACNE-AP, and graph visualization—are currently **pending implementation** in the pipeline. Once implemented, they will be added as new workflows.*

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

- `main.nf` - The primary Nextflow workflow definition containing `ATLAS`, `SURGERY`, `RAYTUNE`, `SUBCLUSTER`, and `PIPELINE`.
- `nextflow.config` - Defines `local`, `local_cpu`, `local_gpu`, and `remote_gpu` executor/container profiles.
- `bin/` - Python scripts executed by Nextflow for the various workflows, including dataset retrieval and model training.
- `configs/` - Base YAML parameter files and the canonical subclustering cohort manifest.
- `containers/subcluster/` - Reproducible FAISS/RAPIDS subclustering image recipe.
- `tests/` - Unit and CPU integration tests for cohort preparation, mapping, and clustering.
- `metadata/` - Contains dataset catalog files like `dsets.csv`.
- `AGENTS.md` and `SKILLS.md` - Context for interacting with the AI agent system to develop or monitor the pipeline.
