# LUCA Single-Cell Coabundance Networks Pipeline

The main objective of this project is to perform an ecological analysis of the cell types in NSCLC lung tumor tissues from scRNA-seq data, in both early and late tumor stages. We obtain coabundance metrics, specifically relying on Mutual Information (MI), as it captures non-linear correlations in the data, and the methodology in the [ARACNE-AP](https://github.com/califano-lab/ARACNe-AP) package ensures a smaller possibility of spurious correlations.

This repository serves as the unified Nextflow automation of the preprocessing, model training, annotation, and downstream analysis. It replaces the old Jupyter notebook workflow with a robust, scalable pipeline.

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

- `main.nf` - The primary Nextflow workflow definition containing the three entry workflows.
- `nextflow.config` - Defines executor profiles (`local`, `remote_gpu`) and singularity bindings.
- `bin/` - Python scripts executed by Nextflow for the various workflows.
- `configs/` - Base YAML configurations for each workflow type.
- `AGENTS.md` and `SKILLS.md` - Context for interacting with the AI agent system to develop or monitor the pipeline.
