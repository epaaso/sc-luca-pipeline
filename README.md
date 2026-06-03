# LUCA Single-Cell Coabundance Networks Pipeline

The main objective of this project is to perform an ecological analysis of the cell types in NSCLC lung tumor tissues from scRNA-seq data, in both early and late tumor stages. We obtain coabundance metrics, specifically relying on Mutual Information (MI), as it captures non-linear correlations in the data, and the methodology in the [ARACNE-AP](https://github.com/califano-lab/ARACNe-AP) package ensures a smaller possibility of spurious correlations.

This repository serves as the unified Nextflow automation of the preprocessing, model training, annotation, and downstream analysis. It replaces the old Jupyter notebook workflow with a robust, scalable pipeline.

![scRNA-seq Workflow](Workflow_dark.png "scRNA-seq Workflow")

## Current Pipeline Status

The pipeline currently implements the deep learning model training and dataset annotation phases:

1. **Atlas Generation (`ATLAS` Workflow)**: Trains the scVI/scANVI reference models to embed and cluster reference datasets.
2. **Dataset Annotation (`SURGERY` Workflow)**: Annotates external datasets via architectural surgery, mapping new cells onto the reference atlas.
3. **Hyperparameter Tuning (`RAYTUNE` Workflow)**: Searches for the best model architecture via Ray Tune.

*Note: The remaining phases—mutual information extraction, network generation with ARACNE-AP, and graph visualization—are currently **pending implementation** in the pipeline. Once implemented, they will be added as new workflows.*

## Execution

This pipeline relies on Nextflow and Singularity for reproducibility. It supports different profiles to easily run on a local workstation or on a remote Slurm GPU cluster.

### Running a Workflow

Use the `-entry` flag to specify which workflow to run: `ATLAS`, `SURGERY`, or `RAYTUNE`. You can override parameters by passing a custom config file with `-c`.

```bash
# Run the dataset surgery workflow using a specific config locally
nextflow run main.nf -entry SURGERY -c configs/surgery_default.yaml -profile local

# Run the atlas generation on a remote Slurm GPU cluster
nextflow run main.nf -entry ATLAS -c configs/atlas_default.yaml -profile remote_gpu
```

### Profiles
- `-profile local`: Runs the processes locally.
- `-profile remote_gpu`: Submits jobs to the Slurm cluster queue.

## Folder Structure

- `main.nf` - The primary Nextflow workflow definition containing the entry workflows (`ATLAS`, `SURGERY`, `RAYTUNE`, and `PIPELINE`).
- `nextflow.config` - Defines executor profiles (`local`, `remote_gpu`) and singularity bindings.
- `bin/` - Python scripts executed by Nextflow for the various workflows, including dataset retrieval and model training.
- `configs/` - Base YAML configurations for each workflow type.
- `metadata/` - Contains dataset catalog files like `dsets.csv`.
- `AGENTS.md` and `SKILLS.md` - Context for interacting with the AI agent system to develop or monitor the pipeline.

## Dataset Catalog (`metadata/dsets.csv`)

The project contains a standardized registry catalog under [metadata/dsets.csv](file:///home/epaaso/REPOS/sc-luca-pipeline/metadata/dsets.csv). This CSV maps and documents the metadata, filtering parameters, and publications for each single-cell RNA-seq dataset integrated into this study.

Key columns in the table include:
- `id`: The unique name identifier of the dataset (e.g., `Chen_Zhang_2020_NSCLC`, `Kim_Lee_2020_LUAD`, `Zuani_2024_NSCLC`).
- `input_adata`: The expected local path to the raw/processed AnnData file.
- `min_counts`, `max_counts`, `min_genes`, `max_genes`, `max_pct_mito`: Quality control (QC) thresholds applied during preprocessing.
- `integrated cells n`: Total cells included from the study.
- `doi`, `journal`, `names`: Citation and publication registry information.
- `stage`, `disease`, `cell_sorting`, `prior_treatment`: Biological details and clinical variables of the dataset.

## Automatic Dataset Downloading & Preprocessing

To make execution seamless and remote-friendly (especially on environments without access to local shared mounts like `/datos`), the pipeline automatically downloads and prepares input files if they do not exist locally:

1. **LUCA Reference Atlas**:
   - If the files `/data/luca_atlas/extended.h5ad` or `/data/luca_atlas/extended_tumor_hvg.h5ad` are missing, Nextflow will trigger `PREPARE_DATASET` for the `extended` or `extended_tumor_hvg` dataset.
   - The pipeline retrieves the full extended atlas directly from CELLxGENE's stable dataset asset portal (`https://datasets.cellxgene.cziscience.com/173984ce-d33a-46b6-ae96-4be47f6c67e8.h5ad`).
   - For `extended_tumor_hvg`, the script automatically subsets the downloaded atlas in memory to `tumor_primary` cells, filters for specific UICC stages (`I` through `IV`), excludes conflicting studies, filters for highly variable genes, and prepares raw counts.

2. **Query/Surgery Datasets**:
   - If the input query file specified for dataset surgery (e.g., `Bishoff`, `Hu`, `Zuani`, `Deng`) does not exist at its designated location, the pipeline automatically fetches the raw matrices/barcodes/features from their public endpoints (EBI BioStudies, Figshare, GEO, CodeOcean) and runs the corresponding filtering pipeline defined in `bin/download_and_preprocess_dataset.py` before executing the surgery mapping.
