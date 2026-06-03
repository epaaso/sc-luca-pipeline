#!/usr/bin/env nextflow

import groovy.json.JsonOutput

nextflow.enable.dsl = 2

def setDefault = { key, value ->
    if (!params.containsKey(key) || params[key] == null) {
        params[key] = value
    }
}

// Global resource defaults
setDefault('slurm_cpus', 8)
setDefault('slurm_memory', '32 GB')
setDefault('slurm_time', '12:00:00')
setDefault('slurm_gpus', 0)

// ---------------------------------------------------------
// ATLAS PROCESS
// ---------------------------------------------------------

process TRAIN_SCANVI_ATLAS {
    label "gpu"
    stageInMode "copy"
    cpus params.slurm_cpus as int
    memory params.slurm_memory
    time params.slurm_time
    clusterOptions { params.slurm_gpus ? "--gres=gpu:rtx5080:${params.slurm_gpus}" : "" }

    input:
    path train_script

    output:
    path "run_config.json"
    path "metrics.json", optional: true
    path "versions.json", optional: true
    path "scvi_history.csv", optional: true
    path "scanvi_history.csv", optional: true
    path "ref_latent.h5ad", optional: true
    path "prepared_atlas.h5ad", optional: true

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    def config = [
        experiment_name: params.experiment_name,
        input_h5ad: params.input_h5ad,
        shared_run_root: params.shared_run_root,
        run_dir: runDir,
        seed: params.seed,
        batch_key: params.batch_key,
        labels_key: params.labels_key,
        unlabeled_category: params.unlabeled_category,
        origin_key: params.origin_key,
        origin_values: params.origin_values,
        stage_key: params.stage_key,
        stages: params.stages,
        study_key: params.study_key,
        exclude_studies: params.exclude_studies,
        sample_key: params.sample_key,
        split_batch_dataset: params.split_batch_dataset,
        counts_layer: params.counts_layer,
        drop_layers: params.drop_layers,
        drop_raw: params.drop_raw,
        round_counts: params.round_counts,
        hvg_mode: params.hvg_mode,
        hvg_n_top_genes: params.hvg_n_top_genes,
        hvg_flavor: params.hvg_flavor,
        hvg_batch_key: params.hvg_batch_key,
        hvg_column: params.hvg_column,
        max_cells: params.max_cells,
        max_genes: params.max_genes,
        save_prepared_h5ad: params.save_prepared_h5ad,
        dl_num_workers: params.dl_num_workers,
        float32_matmul_precision: params.float32_matmul_precision,
        numba_threads: params.numba_threads,
        compute_neighbors: params.compute_neighbors,
        compute_leiden: params.compute_leiden,
        compute_umap: params.compute_umap,
        n_neighbors: params.n_neighbors,
        collapse_wu_dataset: params.collapse_wu_dataset,
        latent_obs_columns: params.latent_obs_columns,
        tensorboard_logging: params.tensorboard_logging,
        tensorboard_log_dir: params.tensorboard_log_dir,
        reuse_scvi_model_if_exists: params.reuse_scvi_model_if_exists,
        scvi_model_params: params.scvi_model_params,
        scvi_train_params: params.scvi_train_params,
        scvi_early_stopping: params.scvi_early_stopping,
        scvi_early_stopping_kwargs: params.scvi_early_stopping_kwargs,
        scvi_plan_kwargs: params.scvi_plan_kwargs,
        scvi_trainer_kwargs: params.scvi_trainer_kwargs,
        scanvi_train_params: params.scanvi_train_params,
        scanvi_early_stopping: params.scanvi_early_stopping,
        scanvi_early_stopping_kwargs: params.scanvi_early_stopping_kwargs,
        scanvi_plan_kwargs: params.scanvi_plan_kwargs,
        scanvi_trainer_kwargs: params.scanvi_trainer_kwargs
    ]
    def configJson = JsonOutput.prettyPrint(JsonOutput.toJson(config))
    """
    set -euo pipefail

    mkdir -p "${runDir}"

    cat > run_config.json <<'JSON'
${configJson}
JSON

    mkdir -p .cache/matplotlib .cache/numba .cache/torch

    singularity exec --nv \\
      --home "\$PWD" \\
      --env XDG_CACHE_HOME="\$PWD/.cache" \\
      --env MPLCONFIGDIR="\$PWD/.cache/matplotlib" \\
      --env NUMBA_CACHE_DIR="\$PWD/.cache/numba" \\
      --env TORCH_HOME="\$PWD/.cache/torch" \\
      --bind "${params.bind_paths}" \\
      "${params.sif}" \\
      python "${train_script}" --config run_config.json

    cp run_config.json "${runDir}/run_config.json"
    for f in metrics.json versions.json scvi_history.csv scanvi_history.csv ref_latent.h5ad prepared_atlas.h5ad; do
      if [[ -f "\$f" ]]; then
        cp "\$f" "${runDir}/\$f"
      fi
    done
    """
}

// ---------------------------------------------------------
// SURGERY PROCESS
// ---------------------------------------------------------

process TRAIN_SCANVI_SURGERY {
    label "gpu"
    stageInMode "copy"
    cpus params.slurm_cpus as int
    memory params.slurm_memory
    time params.slurm_time
    clusterOptions { params.slurm_gpus ? "--gres=gpu:rtx5080:${params.slurm_gpus}" : "" }

    input:
    path train_script

    output:
    path "run_config.json"
    path "query.h5ad", optional: true
    path "query_latent.h5ad", optional: true
    path "*_predicted.csv", optional: true
    path "*_ensembl.csv", optional: true

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    def config = [
        experiment_name: params.experiment_name,
        input_h5ad: params.input_h5ad,
        reference_model: params.reference_model,
        shared_run_root: params.shared_run_root,
        run_dir: runDir,
        seed: params.seed,
        condition_key: params.condition_key,
        cell_type_key: params.cell_type_key,
        dataset_name: params.dataset_name,
        dataset_name_short: params.dataset_name_short,
        dl_num_workers: params.dl_num_workers,
        float32_matmul_precision: params.float32_matmul_precision,
        map_ensembl: params.map_ensembl,
        max_cells: params.max_cells,
        scanvi_train_params: params.scanvi_train_params,
        scanvi_early_stopping: params.scanvi_early_stopping,
        scanvi_early_stopping_kwargs: params.scanvi_early_stopping_kwargs,
        scanvi_plan_kwargs: params.scanvi_plan_kwargs
    ]
    def configJson = JsonOutput.prettyPrint(JsonOutput.toJson(config))
    """
    set -euo pipefail

    mkdir -p "${runDir}"

    cat > run_config.json <<'JSON'
${configJson}
JSON

    mkdir -p .cache/matplotlib .cache/numba .cache/torch

    singularity exec --nv \\
      --home "\$PWD" \\
      --env PYTHONPATH="/home/epaaso/REPOS/sc-luca-explore/utils:\$PYTHONPATH" \\
      --env XDG_CACHE_HOME="\$PWD/.cache" \\
      --env MPLCONFIGDIR="\$PWD/.cache/matplotlib" \\
      --env NUMBA_CACHE_DIR="\$PWD/.cache/numba" \\
      --env TORCH_HOME="\$PWD/.cache/torch" \\
      --bind "${params.bind_paths}" \\
      "${params.sif}" \\
      python "${train_script}" --config run_config.json

    cp run_config.json "${runDir}/run_config.json"
    for f in query.h5ad query_latent.h5ad *_predicted.csv *_ensembl.csv; do
      if [[ -f "\$f" ]]; then
        cp "\$f" "${runDir}/\$f"
      fi
    done
    cp -r surgery_model "${runDir}/" || true
    """
}

// ---------------------------------------------------------
// RAYTUNE PROCESS
// ---------------------------------------------------------

process RUN_SCVI_RAYTUNE {
    label "gpu"
    stageInMode "copy"
    cpus params.slurm_cpus as int
    memory params.slurm_memory
    time params.slurm_time
    clusterOptions { params.slurm_gpus ? "--gres=gpu:rtx5080:${params.slurm_gpus}" : "" }

    input:
    path train_script

    output:
    path "run_config.json"
    path "summary.csv", optional: true
    path "best_config.json", optional: true
    path "best_result.json", optional: true
    path "versions.json", optional: true
    path "ray_log_dir.txt", optional: true

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    def config = [
        experiment_name: params.experiment_name,
        input_h5ad: params.input_h5ad,
        output_root: params.output_root,
        shared_run_root: params.shared_run_root,
        run_dir: runDir,
        metric: params.metric,
        mode: params.mode,
        num_samples: params.num_samples,
        seed: params.seed,
        batch_key: params.batch_key,
        labels_key: params.labels_key,
        model: params.model,
        resources: [
            cpu: params.cpus_per_trial,
            gpu: params.gpus_per_trial
        ],
        ray: [
            num_cpus: params.ray_cpus,
            num_gpus: params.ray_gpus
        ],
        search_space: params.search_space,
        initial_points: params.initial_points,
        scheduler: params.scheduler,
        searcher: params.searcher,
        scheduler_kwargs: params.scheduler_kwargs,
        searcher_kwargs: params.searcher_kwargs,
        save_checkpoints: params.save_checkpoints,
        log_to_driver: params.log_to_driver,
        max_cells: params.max_cells,
        max_genes: params.max_genes
    ]
    def configJson = JsonOutput.prettyPrint(JsonOutput.toJson(config))
    """
    set -euo pipefail

    mkdir -p "${runDir}"

    cat > run_config.json <<'JSON'
${configJson}
JSON

    mkdir -p .cache/matplotlib .cache/numba .cache/ray

    singularity exec --nv \\
      --home "\$PWD" \\
      --env XDG_CACHE_HOME="\$PWD/.cache" \\
      --env MPLCONFIGDIR="\$PWD/.cache/matplotlib" \\
      --env NUMBA_CACHE_DIR="\$PWD/.cache/numba" \\
      --env RAY_TMPDIR="\$PWD/.cache/ray" \\
      --bind "${params.bind_paths}" \\
      "${params.sif}" \\
      python "${train_script}" --config run_config.json

    cp run_config.json "${runDir}/run_config.json"
    for f in summary.csv best_config.json best_result.json versions.json ray_log_dir.txt; do
      if [[ -f "\$f" ]]; then
        cp "\$f" "${runDir}/\$f"
      fi
    done
    """
}

// ---------------------------------------------------------
// WORKFLOW DEFINITIONS
// ---------------------------------------------------------

workflow ATLAS {
    // We only set defaults that are absolutely required if they are missing
    setDefault('experiment_name', 'scanvi_atlas_smoke')
    setDefault('shared_run_root', '/datos/LUCA_model/scanvi_atlas/runs')
    setDefault('sif', '/data/containers/scanvi-atlas-cu12.sif')
    setDefault('bind_paths', '/datos,/data,/home,/datos/home')
    setDefault('train_script', "$baseDir/bin/train_scanvi_atlas.py")
    setDefault('slurm_cpus', 20)
    setDefault('slurm_memory', '96 GB')
    setDefault('slurm_time', '48:00:00')
    setDefault('slurm_gpus', 1)

    new File("${params.shared_run_root}/${params.experiment_name}").mkdirs()
    train_script = file(params.train_script)
    TRAIN_SCANVI_ATLAS(train_script)
}

workflow SURGERY {
    setDefault('experiment_name', 'scanvi_surgery_smoke')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/scanvi_surgery/runs')
    setDefault('sif', '/data/containers/scvi-raytune-py313-cu12.sif')
    setDefault('bind_paths', '/datos,/data,/home,/datos/home')
    setDefault('train_script', "$baseDir/bin/run_surgery.py")
    setDefault('slurm_cpus', 8)
    setDefault('slurm_memory', '64 GB')
    setDefault('slurm_time', '12:00:00')
    setDefault('slurm_gpus', 1)

    new File("${params.shared_run_root}/${params.experiment_name}").mkdirs()
    train_script = file(params.train_script)
    TRAIN_SCANVI_SURGERY(train_script)
}

workflow RAYTUNE {
    setDefault('experiment_name', 'scvi_raytune_smoke')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/raytune_runs')
    setDefault('sif', '/data/containers/scvi-raytune-cu12.sif')
    setDefault('bind_paths', '/datos,/data,/home,/datos/home')
    setDefault('train_script', "$baseDir/bin/train_scvi_raytune.py")
    setDefault('slurm_cpus', 4)
    setDefault('slurm_memory', '32 GB')
    setDefault('slurm_time', '02:00:00')
    setDefault('slurm_gpus', 1)
    
    // Defaulting these params needed for raytune config assembly
    setDefault('metric', 'validation_loss')
    setDefault('mode', 'min')
    setDefault('num_samples', 1)
    setDefault('seed', 0)
    setDefault('batch_key', 'dataset')
    setDefault('labels_key', 'cell_type_tumor')
    setDefault('model', 'SCVI')
    setDefault('cpus_per_trial', 4)
    setDefault('gpus_per_trial', 1)
    setDefault('ray_cpus', 4)
    setDefault('ray_gpus', 1)

    new File("${params.shared_run_root}/${params.experiment_name}").mkdirs()
    train_script = file(params.train_script)
    RUN_SCVI_RAYTUNE(train_script)
}
