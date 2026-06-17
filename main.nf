#!/usr/bin/env nextflow

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

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
setDefault('subcluster_backend', 'cpu')
setDefault('subcluster_container', 'sc-luca-subcluster:latest')
setDefault('subcluster_cpus', 16)
setDefault('subcluster_memory', '64 GB')
setDefault('subcluster_time', '12:00:00')
setDefault('subcluster_gpus', 0)
setDefault('subcluster_python', 'python')
setDefault('subcluster_pythonpath', '')
setDefault('subcluster_ld_library_path', '')
setDefault('graph_container', 'sc-luca-subcluster:latest')
setDefault('graph_tool_container', 'tiagopeixoto/graph-tool:latest')
setDefault('graph_python', 'python')
setDefault('graph_pythonpath', '')
setDefault('graph_ld_library_path', '')
setDefault('graph_cpus', 8)
setDefault('graph_memory', '32 GB')
setDefault('graph_time', '24:00:00')

def subclusterDefaults = {
    setDefault('subcluster_script', "$baseDir/bin/run_subcluster.py")
    setDefault('subcluster_input_h5ad', null)
    setDefault(
        'subcluster_cohort_manifest',
        params.subcluster_input_h5ad ? null : "$baseDir/configs/subcluster_cohort_default.json"
    )
    setDefault('subcluster_engine', 'faiss')
    setDefault('subcluster_backend', 'cpu')
    setDefault('subcluster_gpu_postprocess', 'cpu')
    setDefault('subcluster_container', 'sc-luca-subcluster:latest')
    setDefault('subcluster_seed', 42)
    setDefault('subcluster_n_neighbors', 30)
    setDefault('subcluster_leiden_resolution', 0.5)
    setDefault('subcluster_umap_min_dist', 0.03)
    setDefault('subcluster_umap_spread', 0.6)
    setDefault('subcluster_faiss_nlist', null)
    setDefault('subcluster_faiss_nprobe', null)
    setDefault('subcluster_tumor_pattern', 'Tumor')
    setDefault('subcluster_epithelial_types', [
        'Alveolar cell type 1',
        'Alveolar cell type 2',
        'ROS1+ healthy epithelial',
        'transitional club/AT2',
        'Club',
        'Ciliated'
    ])
    setDefault('early_mapping_json', null)
    setDefault('late_mapping_json', null)
    setDefault('subcluster_cpus', 16)
    setDefault('subcluster_memory', '64 GB')
    setDefault('subcluster_time', '12:00:00')
    setDefault('subcluster_gpus', 0)
    setDefault('subcluster_python', 'python')
    setDefault('subcluster_pythonpath', '')
    setDefault('subcluster_ld_library_path', '')
}

def graphDefaults = {
    setDefault('graph_script', "$baseDir/bin/run_graph_phase.py")
    setDefault('graph_input', null)
    setDefault('graph_step', 'all')
    setDefault('graph_global_dir', null)
    setDefault('ecotype_dir', null)
    setDefault('graph_cell_type_key', 'cell_type_adjusted')
    setDefault('graph_sample_key', 'batch')
    setDefault('graph_dataset_key', 'dataset')
    setDefault('graph_stage_key', 'stage')
    setDefault('graph_min_cells_per_type_sample', 2)
    setDefault('graph_seed', 42)
    setDefault('ecotype_k', 4)
    setDefault('ecotype_resolution', 0.2)
    setDefault('aracne_jar', '/data/containers/ARACNe-AP/dist/aracne.jar')
    setDefault('aracne_java', '/data/containers/java/bin/java')
    setDefault('aracne_pvalue', '1E-8')
    setDefault('aracne_bootstraps', 500)
    setDefault('aracne_java_memory', '8G')
    setDefault('aracne_min_samples', 5)
    setDefault('aracne_min_nonzero_types', 3)
    setDefault('aracne_dry_run', false)
}

def absoluteParamPath = { value ->
    value ? file(value).toAbsolutePath().toString() : null
}

def subclusterConfig = { manifestPath, inputPath ->
    [
        subcluster_cohort_manifest: manifestPath,
        subcluster_input_h5ad: inputPath,
        subcluster_output_dir: 'subcluster_results',
        subcluster_engine: params.subcluster_engine,
        subcluster_backend: params.subcluster_backend,
        subcluster_gpu_postprocess: params.subcluster_gpu_postprocess,
        subcluster_seed: params.subcluster_seed,
        subcluster_n_neighbors: params.subcluster_n_neighbors,
        subcluster_leiden_resolution: params.subcluster_leiden_resolution,
        subcluster_umap_min_dist: params.subcluster_umap_min_dist,
        subcluster_umap_spread: params.subcluster_umap_spread,
        subcluster_faiss_nlist: params.subcluster_faiss_nlist,
        subcluster_faiss_nprobe: params.subcluster_faiss_nprobe,
        subcluster_tumor_pattern: params.subcluster_tumor_pattern,
        subcluster_epithelial_types: params.subcluster_epithelial_types
    ]
}

// ---------------------------------------------------------
// DATA PREPARATION / DOWNLOAD
// ---------------------------------------------------------

process PREPARE_DATASET {
    stageInMode "copy"
    cpus 4
    memory '16 GB'

    input:
    val dataset_name
    path python_script

    output:
    path "prepared_dataset.h5ad"

    script:
    """
    python "${python_script}" --dataset "${dataset_name}" --output "prepared_dataset.h5ad"
    """
}

// ---------------------------------------------------------
// CONFIGURATION GENERATORS
// ---------------------------------------------------------

process PREPARE_ATLAS_CONFIG {
    executor 'local'

    output:
    path "run_config.json"

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    def config = [
        experiment_name: params.experiment_name,
        input_h5ad: "atlas_input.h5ad",
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
    cat > run_config.json <<'JSON'
${configJson}
JSON
    """
}

process PREPARE_SURGERY_CONFIG {
    executor 'local'

    input:
    val reference_model_path

    output:
    path "run_config.json"

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    def config = [
        experiment_name: params.experiment_name,
        input_h5ad: "query.h5ad",
        reference_model: reference_model_path,
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
    cat > run_config.json <<'JSON'
${configJson}
JSON
    """
}

process APPLY_RAYTUNE_CONFIG {
    executor 'local'

    input:
    path base_config
    path best_config
    path apply_script

    output:
    path "run_config.yaml"

    script:
    """
    python "${apply_script}" --base "${base_config}" --best-config "${best_config}" --output run_config.yaml
    """
}

process PREPARE_SUBCLUSTER_CONFIG {
    executor 'local'

    input:
    path subcluster_input

    output:
    path "subcluster_config.json", emit: config_file

    script:
    def manifestPath = params.subcluster_input_h5ad ? null : subcluster_input.toString()
    def inputPath = params.subcluster_input_h5ad ? subcluster_input.toString() : null
    def configJson = JsonOutput.prettyPrint(JsonOutput.toJson(subclusterConfig(manifestPath, inputPath)))
    """
    cat > subcluster_config.json <<'JSON'
${configJson}
JSON
    """
}

process PREPARE_SUBCLUSTER_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path subcluster_script

    output:
    path "run_subcluster.py", emit: script_file

    script:
    """
    test -f "${subcluster_script}"
    """
}

process PREPARE_GRAPH_CONFIG {
    executor 'local'

    input:
    path graph_input

    output:
    path "graph_config.json"

    script:
    def config = [
        graph_input: graph_input.toString(),
        graph_output_dir: 'graph_results',
        graph_global_dir: params.graph_global_dir,
        ecotype_dir: params.ecotype_dir,
        graph_cell_type_key: params.graph_cell_type_key,
        graph_sample_key: params.graph_sample_key,
        graph_dataset_key: params.graph_dataset_key,
        graph_stage_key: params.graph_stage_key,
        graph_min_cells_per_type_sample: params.graph_min_cells_per_type_sample,
        graph_seed: params.graph_seed,
        ecotype_k: params.ecotype_k,
        ecotype_resolution: params.ecotype_resolution,
        aracne_jar: params.aracne_jar,
        aracne_java: params.aracne_java,
        aracne_pvalue: params.aracne_pvalue,
        aracne_bootstraps: params.aracne_bootstraps,
        aracne_java_memory: params.aracne_java_memory,
        aracne_min_samples: params.aracne_min_samples,
        aracne_min_nonzero_types: params.aracne_min_nonzero_types,
        aracne_dry_run: params.aracne_dry_run
    ]
    def configJson = JsonOutput.prettyPrint(JsonOutput.toJson(config))
    """
    cat > graph_config.json <<'JSON'
${configJson}
JSON
    """
}

process PREPARE_GRAPH_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path graph_script

    output:
    path "run_graph_phase.py"

    script:
    """
    test -f "${graph_script}"
    """
}

process PREPARE_COHORT_QUERY {
    executor 'local'
    stageInMode "copy"
    cpus 4
    memory '16 GB'

    input:
    val entry
    path download_script

    output:
    tuple val(entry), path("${entry.name}_input.h5ad")

    script:
    def sourcePath = entry.input_h5ad
    def downloadName = entry.download_name ?: entry.name
    """
    set -euo pipefail
    if [[ -f "${sourcePath}" ]]; then
      cp "${sourcePath}" "${entry.name}_input.h5ad"
    else
      python "${download_script}" --dataset "${downloadName}" --output "${entry.name}_input.h5ad"
    fi
    """
}

process PREPARE_PIPELINE_SUBCLUSTER_CONFIG {
    executor 'local'

    input:
    path atlas_latent
    path cohort_dirs

    output:
    path "pipeline_cohort.json", emit: cohort_manifest
    path "subcluster_config.json", emit: config_file

    script:
    def configJson = JsonOutput.prettyPrint(JsonOutput.toJson(subclusterConfig('pipeline_cohort.json', null)))
    def cohortArgs = cohort_dirs.collect { "\"${it}\"" }.join(' ')
    """
    set -euo pipefail
    python - "${atlas_latent}" ${cohortArgs} <<'PY'
import json
import pathlib
import sys

atlas = pathlib.Path(sys.argv[1]).resolve()
queries = []
for directory in sys.argv[2:]:
    directory = pathlib.Path(directory).resolve()
    with (directory / "entry.json").open() as handle:
        entry = json.load(handle)
    entry["latent_h5ad"] = str(directory / "query_latent.h5ad")
    queries.append(entry)

manifest = {
    "atlas": {
        "name": "Atlas",
        "latent_h5ad": str(atlas),
        "representation": "X",
        "cell_type_key": "cell_type",
        "batch_key": "batch",
        "stage_key": "uicc_stage"
    },
    "queries": queries
}
with open("pipeline_cohort.json", "w") as handle:
    json.dump(manifest, handle, indent=2)
PY
    cat > subcluster_config.json <<'JSON'
${configJson}
JSON
    """
}

// ---------------------------------------------------------
// ML TRAINING PROCESSES
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
    path 'raytune_input.h5ad'

    output:
    path "run_config.json"
    path "summary.csv", optional: true
    path "best_config.json", optional: true, emit: best_config_json
    path "best_result.json", optional: true
    path "versions.json", optional: true
    path "ray_log_dir.txt", optional: true

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    def config = [
        experiment_name: params.experiment_name,
        input_h5ad: "raytune_input.h5ad",
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

process TRAIN_SCANVI_ATLAS {
    label "gpu"
    stageInMode "copy"
    cpus params.slurm_cpus as int
    memory params.slurm_memory
    time params.slurm_time
    clusterOptions { params.slurm_gpus ? "--gres=gpu:rtx5080:${params.slurm_gpus}" : "" }

    input:
    path train_script
    path config_file
    path 'atlas_input.h5ad'

    output:
    path "${config_file}"
    path "metrics.json", optional: true
    path "versions.json", optional: true
    path "scvi_history.csv", optional: true
    path "scanvi_history.csv", optional: true
    path "ref_latent.h5ad", optional: true, emit: ref_latent_h5ad
    path "prepared_atlas.h5ad", optional: true
    path "scanvi_model", optional: true, emit: scanvi_model

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    """
    set -euo pipefail

    mkdir -p "${runDir}"
    mkdir -p .cache/matplotlib .cache/numba .cache/torch

    singularity exec --nv \\
      --home "\$PWD" \\
      --env XDG_CACHE_HOME="\$PWD/.cache" \\
      --env MPLCONFIGDIR="\$PWD/.cache/matplotlib" \\
      --env NUMBA_CACHE_DIR="\$PWD/.cache/numba" \\
      --env TORCH_HOME="\$PWD/.cache/torch" \\
      --bind "${params.bind_paths}" \\
      "${params.sif}" \\
      python "${train_script}" --config "${config_file}"

    cp "${config_file}" "${runDir}/run_config.json"
    for f in metrics.json versions.json scvi_history.csv scanvi_history.csv ref_latent.h5ad prepared_atlas.h5ad; do
      if [[ -f "\$f" ]]; then
        cp "\$f" "${runDir}/\$f"
      fi
    done

    if [[ -d "${runDir}/scanvi_model" ]]; then
      cp -r "${runDir}/scanvi_model" ./scanvi_model
    fi
    """
}

process TRAIN_SCANVI_SURGERY {
    label "gpu"
    stageInMode "copy"
    cpus params.slurm_cpus as int
    memory params.slurm_memory
    time params.slurm_time
    clusterOptions { params.slurm_gpus ? "--gres=gpu:rtx5080:${params.slurm_gpus}" : "" }

    input:
    path train_script
    path config_file
    path 'query.h5ad'
    path 'scanvi_model'

    output:
    path "${config_file}"
    path "query.h5ad", optional: true
    path "query_latent.h5ad", optional: true
    path "*_predicted.csv", optional: true
    path "*_ensembl.csv", optional: true

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}"
    """
    set -euo pipefail

    mkdir -p "${runDir}"
    mkdir -p .cache/matplotlib .cache/numba .cache/torch

    singularity exec --nv \\
      --home "\$PWD" \\
      --env PYTHONPATH="/home/epaaso/REPOS/sc-luca-explore/utils:\${PYTHONPATH:-}" \\
      --env XDG_CACHE_HOME="\$PWD/.cache" \\
      --env MPLCONFIGDIR="\$PWD/.cache/matplotlib" \\
      --env NUMBA_CACHE_DIR="\$PWD/.cache/numba" \\
      --env TORCH_HOME="\$PWD/.cache/torch" \\
      --bind "${params.bind_paths}" \\
      "${params.sif}" \\
      python "${train_script}" --config "${config_file}"

    cp "${config_file}" "${runDir}/run_config.json"
    for f in query.h5ad query_latent.h5ad *_predicted.csv *_ensembl.csv; do
      if [[ -f "\$f" ]]; then
        cp "\$f" "${runDir}/\$f"
      fi
    done
    if [[ -d scanvi_model ]]; then
      cp -r scanvi_model "${runDir}/"
    fi
    """
}

process RUN_COHORT_SURGERY {
    label "gpu"
    stageInMode "copy"
    cpus params.slurm_cpus as int
    memory params.slurm_memory
    time params.slurm_time
    clusterOptions { params.slurm_gpus ? "--gres=gpu:rtx5080:${params.slurm_gpus}" : "" }

    input:
    tuple val(entry), path(query_file)
    path 'scanvi_model'
    path surgery_script

    output:
    path "${entry.name}", emit: cohort_dir

    script:
    def name = entry.name as String
    def runDir = "${params.shared_run_root}/${params.experiment_name}/surgery/${name}"
    def config = [
        experiment_name: "${params.experiment_name}_${name}",
        input_h5ad: query_file.toString(),
        reference_model: "scanvi_model",
        shared_run_root: params.shared_run_root,
        run_dir: runDir,
        seed: params.seed,
        condition_key: params.condition_key,
        cell_type_key: params.cell_type_key,
        dataset_name: name,
        dataset_name_short: name,
        dl_num_workers: params.dl_num_workers,
        float32_matmul_precision: params.float32_matmul_precision,
        map_ensembl: params.map_ensembl,
        max_cells: params.max_cells,
        scanvi_train_params: params.scanvi_train_params,
        scanvi_early_stopping: params.scanvi_early_stopping,
        scanvi_early_stopping_kwargs: params.scanvi_early_stopping_kwargs,
        scanvi_plan_kwargs: params.scanvi_plan_kwargs
    ]
    def entryOutput = [
        name: name,
        representation: 'X_scVI',
        cell_type_key: 'predicted_cell_type',
        batch_key: entry.batch_key ?: 'sample',
        stage_key: entry.stage_key,
        stage_override: entry.stage_override,
        stage_replacements: entry.stage_replacements
    ].findAll { key, value -> value != null }
    def configJson = JsonOutput.prettyPrint(JsonOutput.toJson(config))
    def entryJson = JsonOutput.prettyPrint(JsonOutput.toJson(entryOutput))
    """
    set -euo pipefail

    mkdir -p "${runDir}" "${name}" .cache/matplotlib .cache/numba .cache/torch
    cat > run_config.json <<'JSON'
${configJson}
JSON

    singularity exec --nv \\
      --home "\$PWD" \\
      --env XDG_CACHE_HOME="\$PWD/.cache" \\
      --env MPLCONFIGDIR="\$PWD/.cache/matplotlib" \\
      --env NUMBA_CACHE_DIR="\$PWD/.cache/numba" \\
      --env TORCH_HOME="\$PWD/.cache/torch" \\
      --bind "${params.bind_paths}" \\
      "${params.sif}" \\
      python "${surgery_script}" --config run_config.json

    cp "${runDir}/query_latent.h5ad" "${name}/query_latent.h5ad"
    cat > "${name}/entry.json" <<'JSON'
${entryJson}
JSON
    """
}

process RUN_SUBCLUSTER {
    label "subcluster"
    container params.subcluster_container
    stageInMode "copy"
    cpus params.subcluster_cpus as int
    memory params.subcluster_memory
    time params.subcluster_time
    clusterOptions { params.subcluster_backend == 'gpu' && params.subcluster_gpus ? "--gres=gpu:rtx5080:${params.subcluster_gpus}" : "" }

    input:
    path subcluster_script
    path config_file
    path cohort_manifest

    output:
    path "subcluster_results", emit: raw_results

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}/subcluster"
    """
    set -euo pipefail
    mkdir -p .cache/numba .cache/matplotlib
    export NUMBA_NUM_THREADS="${params.subcluster_cpus}"
    export NUMBA_CACHE_DIR="\$PWD/.cache/numba"
    export MPLCONFIGDIR="\$PWD/.cache/matplotlib"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export PYTHONPATH="${params.subcluster_pythonpath}${params.subcluster_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.subcluster_ld_library_path}${params.subcluster_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    ${params.subcluster_backend == 'gpu' ? 'nvidia-smi -L' : 'true'}
    "${params.subcluster_python}" "${subcluster_script}" --config "${config_file}"
    mkdir -p "${runDir}"
    cp -rn subcluster_results/. "${runDir}/"
    """

    stub:
    """
    mkdir -p subcluster_results
    touch subcluster_results/stub.txt
    """
}

process APPLY_SUBCLUSTER_MAPPINGS {
    label "subcluster"
    container params.subcluster_container
    stageInMode "copy"
    cpus 2
    memory '8 GB'
    time '02:00:00'

    input:
    path subcluster_script
    path raw_results
    val early_mapping_path
    val late_mapping_path

    output:
    path "mapping_results"

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}/subcluster"
    def mappingConfigJson = JsonOutput.prettyPrint(JsonOutput.toJson([
        early_mapping_json: early_mapping_path ?: null,
        late_mapping_json: late_mapping_path ?: null
    ]))
    """
    set -euo pipefail
    export PYTHONPATH="${params.subcluster_pythonpath}${params.subcluster_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.subcluster_ld_library_path}${params.subcluster_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    cat > subcluster_mapping_config.json <<'JSON'
${mappingConfigJson}
JSON
    "${params.subcluster_python}" "${subcluster_script}" \\
      --config subcluster_mapping_config.json \\
      --apply-mappings-from "${raw_results}" \\
      --mapping-output-dir mapping_results
    mkdir -p "${runDir}"
    cp -rn mapping_results/. "${runDir}/"
    """

    stub:
    """
    mkdir -p mapping_results
    touch mapping_results/mapping_status.json
    """
}

process RUN_GRAPH_PHASE {
    label "graph"
    container params.graph_container
    stageInMode "copy"
    cpus params.graph_cpus as int
    memory params.graph_memory
    time params.graph_time

    input:
    path graph_script
    path graph_config
    path graph_input

    output:
    path "graph_results", emit: results

    script:
    def runDir = "${params.shared_run_root}/${params.experiment_name}/graph"
    """
    set -euo pipefail
    export PYTHONPATH="${params.graph_pythonpath}${params.graph_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.graph_ld_library_path}${params.graph_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    "${params.graph_python}" "${graph_script}" --config "${graph_config}" --step "${params.graph_step}"
    mkdir -p "${runDir}"
    cp -rn graph_results/. "${runDir}/"
    """

    stub:
    """
    mkdir -p graph_results
    touch graph_results/stub.txt
    """
}

process RUN_MARKER_EXTRACTION {
    label "graph"
    container params.graph_container
    stageInMode "copy"
    cpus 6
    memory '32 GB'
    time '24:00:00'

    input:
    path marker_script
    path graph_results
    val adata_dir
    val stage
    path modal_de_script

    output:
    tuple val(stage), path("marker_results_${stage}"), emit: results

    script:
    def membership_csv = "${graph_results}/ecotype/membership_${stage}.csv"
    """
    set -euo pipefail
    export PYTHONPATH="${params.graph_pythonpath}${params.graph_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.graph_ld_library_path}${params.graph_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    export NUMBA_CACHE_DIR="\$PWD/.cache/numba"
    "${params.graph_python}" "${marker_script}" \\
      --config "${graph_results}/run_config.json" \\
      --adata-dir "${adata_dir}" \\
      --membership "${membership_csv}" \\
      --stage "${stage}" \\
      --output_dir "marker_results_${stage}"
    """
}

process PREPARE_MARKER_ANALYSIS_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path script_file

    output:
    path "run_marker_analysis.py"

    script:
    """
    test -f "${script_file}"
    """
}

process PREPARE_GRAPH_ANALYSIS_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path script_file

    output:
    path "run_graph_analysis.py"

    script:
    """
    test -f "${script_file}"
    """
}

process PREPARE_GSEA_GMT {
    executor 'local'
    stageInMode "copy"

    input:
    path gmt_file

    output:
    path "h.all.v2023.2.Hs.symbols.gmt"

    script:
    """
    test -f "${gmt_file}"
    """
}

process PREPARE_CELLPHONEDB_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path script_file

    output:
    path "run_cellphonedb.py"

    script:
    """
    test -f "${script_file}"
    """
}

process PREPARE_CIRCOS_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path script_file

    output:
    path "run_circos_plot.py"

    script:
    """
    test -f "${script_file}"
    """
}

process PREPARE_MARKER_EXTRACTION_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path script_file

    output:
    path "run_marker_extraction.py"

    script:
    """
    test -f "${script_file}"
    """
}

process PREPARE_CELL_MAPPINGS {
    executor 'local'
    stageInMode "copy"

    input:
    path json_file

    output:
    path "cell_mappings.json"

    script:
    """
    test -f "${json_file}"
    """
}

process PREPARE_MODAL_DE_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path script_file

    output:
    path "modal_DE.py"

    script:
    """
    test -f "${script_file}"
    """
}

process PREPARE_CONSENSUS_SCRIPT {
    executor 'local'
    stageInMode "copy"

    input:
    path script_file

    output:
    path "aggregate_consensus.py"

    script:
    """
    chmod +x "${script_file}"
    """
}

process RUN_CONSENSUS_AGGREGATION {
    label "graph"
    container params.graph_container
    stageInMode "copy"
    cpus 1
    memory '8 GB'

    input:
    path aggregate_script
    path marker_root
    val plots_dirs

    output:
    path "consensus_output", emit: results

    script:
    def dirs_string = plots_dirs.join(' ')
    """
    set -euo pipefail
    mkdir -p consensus_output
    export PYTHONPATH="${baseDir}/bin:/datos/home/epaaso/slurm-gpu-jobs/python_packages:${params.graph_pythonpath}${params.graph_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.graph_ld_library_path}${params.graph_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    
    mkdir -p unified_plots/gseapy_gsea
    for d in ${dirs_string}; do
      if [[ -d "\${d}" ]]; then
        cp -rn "\${d}"/. unified_plots/ || true
      fi
    done

    "${params.graph_python}" "${aggregate_script}" \
      --marker-dir "${marker_root}" \
      --gsea-dir "unified_plots/gseapy_gsea" \
      --output-dir "consensus_output"
    """
}


process RUN_MARKER_ANALYSIS {
    label "graph"
    container params.graph_container
    stageInMode "copy"
    cpus 1
    memory '8 GB'

    input:
    path marker_analysis_script
    path marker_file
    path gsea_gmt
    path modal_de_script

    output:
    path "marker_analysis_output", emit: results

    script:
    def name = marker_file.baseName
    def time = name.contains("I-II") ? "I-II" : "III-IV"
    """
    set -euo pipefail
    export PYTHONPATH="${baseDir}/bin:/datos/home/epaaso/slurm-gpu-jobs/python_packages:${params.graph_pythonpath}${params.graph_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.graph_ld_library_path}${params.graph_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    export NUMBA_CACHE_DIR="\$PWD/.cache/numba"
    mkdir -p marker_analysis_output
    "${params.graph_python}" "${marker_analysis_script}" \
      --marker-file "${marker_file}" \
      --gsea-gmt "${gsea_gmt}" \
      --output-dir "marker_analysis_output" \
      --name "${name}" \
      --time "${time}"
    """
}

process RUN_CELLPHONEDB {
    label "graph"
    container params.graph_container
    stageInMode "copy"
    cpus 7
    memory '96 GB'
    time '12:00:00'

    input:
    path cpdb_script
    path "early"
    path "late"
    path graph_results
    val adata_dir
    val cpdb_file
    path modal_de_script

    output:
    path "cpdb_results", emit: results

    script:
    def membership_early = "${graph_results}/ecotype/membership_early.csv"
    def membership_late = "${graph_results}/ecotype/membership_late.csv"
    """
    set -euo pipefail
    echo "Rerunning CellPhoneDB with memory optimization and numpy.rec shim"
    export PYTHONPATH="${baseDir}/bin:/datos/home/epaaso/slurm-gpu-jobs/python_packages:${params.graph_pythonpath}${params.graph_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.graph_ld_library_path}${params.graph_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    export NUMBA_CACHE_DIR="\$PWD/.cache/numba"
    
    mkdir -p merged_markers
    cp -rn early/. merged_markers/ || true
    cp -rn late/. merged_markers/ || true

    "${params.graph_python}" "${cpdb_script}" all \\
      --marker-root "merged_markers" \\
      --output-root "cpdb_results" \\
      --data-dir "${adata_dir}" \\
      --graph-dir "${graph_results}/ecotype_graphs" \\
      --cpdb-file "${cpdb_file}" \\
      --membership-early "${membership_early}" \\
      --membership-late "${membership_late}"
    """
}

process RUN_GRAPH_ANALYSIS {
    label "graph"
    container params.graph_container
    stageInMode "copy"
    cpus 1
    memory '8 GB'

    input:
    path analysis_script
    path graph_results

    output:
    path "analysis_output", emit: results

    script:
    """
    set -euo pipefail
    export PYTHONPATH="${baseDir}/bin:/datos/home/epaaso/slurm-gpu-jobs/python_packages:${params.graph_pythonpath}${params.graph_pythonpath ? ':' : ''}\${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${params.graph_ld_library_path}${params.graph_ld_library_path ? ':' : ''}\${LD_LIBRARY_PATH:-}"
    mkdir -p analysis_output
    
    for stage in early late; do
      net=\$(find "${graph_results}/global" -name "net_\${stage}_MI_pearson.txt" | head -n 1)
      if [[ -n "\${net}" && -f "\${net}" ]]; then
        "${params.graph_python}" "${analysis_script}" \
          --graph-file "\${net}" \
          --output-dir "analysis_output/global/\${stage}" \
          --time-label "\${stage}" \
          --abundance-file "${graph_results}/ecotype/abundance_\${stage}.csv" \
          --global-membership-file "${graph_results}/ecotype/membership_\${stage}.csv"
      fi
    done
    
    for stage in early late; do
      if [[ -d "${graph_results}/ecotype_graphs/\${stage}" ]]; then
        for cluster_dir in "${graph_results}/ecotype_graphs/\${stage}"/cluster_*; do
          if [[ -d "\${cluster_dir}" ]]; then
            cluster=\$(basename "\${cluster_dir}")
            net=\$(find "\${cluster_dir}" -name "*_MI_pearson.txt" | head -n 1)
            if [[ -n "\${net}" && -f "\${net}" ]]; then
              "${params.graph_python}" "${analysis_script}" \
                --graph-file "\${net}" \
                --output-dir "analysis_output/ecotype_graphs/\${stage}/\${cluster}" \
                --time-label "\${stage}_\${cluster}" \
                --groups-file "${graph_results}/global/groups_\${stage}.csv" \
                --membership-file "${graph_results}/ecotype/membership_\${stage}.csv"
            fi
          fi
        done
      fi
    done
    """
}

process RUN_CIRCOS_PLOT {
    label "graph"
    container params.graph_tool_container
    stageInMode "copy"
    cpus 1
    memory '8 GB'

    input:
    path circos_script
    path graph_results
    path metadata

    output:
    path "circos_output", emit: results

    script:
    """
    set -euo pipefail
    mkdir -p circos_output
    
    for stage in early late; do
      net=\$(find "${graph_results}/global" -name "net_\${stage}_MI_pearson.txt" | head -n 1)
      if [[ -n "\${net}" && -f "\${net}" ]]; then
        python "${circos_script}" --graph-file "\${net}" --metadata "${metadata}" --output-dir "circos_output/global/\${stage}" --time-label "\${stage}"
      fi
    done
    
    for stage in early late; do
      if [[ -d "${graph_results}/ecotype_graphs/\${stage}" ]]; then
        for cluster_dir in "${graph_results}/ecotype_graphs/\${stage}"/cluster_*; do
          if [[ -d "\${cluster_dir}" ]]; then
            cluster=\$(basename "\${cluster_dir}")
            net=\$(find "\${cluster_dir}" -name "*_MI_pearson.txt" | head -n 1)
            if [[ -n "\${net}" && -f "\${net}" ]]; then
              python "${circos_script}" --graph-file "\${net}" --metadata "${metadata}" --output-dir "circos_output/ecotype_graphs/\${stage}/\${cluster}" --time-label "\${stage}_\${cluster}"
            fi
          fi
        done
      fi
    done
    """
}

// ---------------------------------------------------------
// WORKFLOW DEFINITIONS
// ---------------------------------------------------------

workflow ATLAS {
    setDefault('experiment_name', 'scanvi_atlas_smoke')
    setDefault('shared_run_root', '/datos/LUCA_model/scanvi_atlas/runs')
    setDefault('sif', '/data/containers/scanvi-atlas-cu12.sif')
    setDefault('bind_paths', '/datos,/data,/home,/datos/home')
    setDefault('train_script', "$baseDir/bin/train_scanvi_atlas.py")
    setDefault('download_script', "$baseDir/bin/download_and_preprocess_dataset.py")
    setDefault('input_h5ad', '/data/luca_atlas/extended.h5ad')
    setDefault('slurm_cpus', 20)
    setDefault('slurm_memory', '96 GB')
    setDefault('slurm_time', '48:00:00')
    setDefault('slurm_gpus', 1)

    // ML defaults
    setDefault('seed', 0)
    setDefault('batch_key', 'dataset')
    setDefault('labels_key', 'cell_type_tumor')
    setDefault('unlabeled_category', 'Unknown')
    setDefault('origin_key', 'origin')
    setDefault('origin_values', ['tumor_primary'])
    setDefault('stage_key', 'uicc_stage')
    setDefault('stages', ['I', 'II', 'III', 'III or IV', 'IV'])
    setDefault('study_key', 'study')
    setDefault('exclude_studies', ['Goveia_Carmeliet_2020', 'Leader_Merad_2021', 'Guo_Zhang_2018', 'Wu_Zhou_2021'])
    setDefault('sample_key', 'sample')
    setDefault('split_batch_dataset', null)
    setDefault('counts_layer', 'count')
    setDefault('drop_layers', ['counts_length_scaled', 'count'])
    setDefault('drop_raw', true)
    setDefault('round_counts', false)
    setDefault('hvg_mode', 'compute')
    setDefault('hvg_n_top_genes', 6000)
    setDefault('hvg_flavor', 'seurat_v3')
    setDefault('hvg_batch_key', 'dataset')
    setDefault('hvg_column', 'is_highly_variable')
    setDefault('max_cells', null)
    setDefault('max_genes', null)
    setDefault('save_prepared_h5ad', false)
    setDefault('dl_num_workers', 8)
    setDefault('float32_matmul_precision', 'high')
    setDefault('numba_threads', 30)
    setDefault('compute_neighbors', true)
    setDefault('compute_leiden', true)
    setDefault('compute_umap', true)
    setDefault('n_neighbors', 8)
    setDefault('collapse_wu_dataset', true)
    setDefault('latent_obs_columns', ['dataset', 'study', 'sample', 'uicc_stage'])
    setDefault('tensorboard_logging', true)
    setDefault('tensorboard_log_dir', null)
    setDefault('reuse_scvi_model_if_exists', true)
    setDefault('scvi_model_params', [n_layers: 4, n_latent: 10, n_hidden: 1024, gene_likelihood: 'nb', dispersion: 'gene-batch', use_batch_norm: 'both', encode_covariates: true, deeply_inject_covariates: false])
    setDefault('scvi_train_params', [max_epochs: 100, batch_size: 128])
    setDefault('scvi_early_stopping', true)
    setDefault('scvi_early_stopping_kwargs', [early_stopping_monitor: 'elbo_validation', early_stopping_patience: 10, early_stopping_min_delta: 0.1])
    setDefault('scvi_plan_kwargs', [reduce_lr_on_plateau: true, lr_patience: 8, lr_factor: 0.05])
    setDefault('scvi_trainer_kwargs', [gradient_clip_val: 10.0])
    setDefault('scanvi_train_params', [max_epochs: 300, batch_size: 128])
    setDefault('scanvi_early_stopping', true)
    setDefault('scanvi_early_stopping_kwargs', [early_stopping_monitor: 'elbo_validation', early_stopping_patience: 10, early_stopping_min_delta: 0.7])
    setDefault('scanvi_plan_kwargs', [reduce_lr_on_plateau: true, lr_patience: 8, lr_factor: 0.1, lr: 1e-4])
    setDefault('scanvi_trainer_kwargs', [gradient_clip_val: 10.0])

    new File("${params.shared_run_root}/${params.experiment_name}").mkdirs()
    PREPARE_ATLAS_CONFIG()
    
    // Determine atlas input dataset (download if not exists)
    def input_file_path = params.input_h5ad
    def atlas_input
    def is_atlas_dataset = (input_file_path =~ /(?i)(extended_tumor_hvg|extended)/)
    def file_exists = file(input_file_path).exists()
    
    if (!file_exists && is_atlas_dataset) {
        def ds_name = ""
        if (input_file_path =~ /(?i)extended_tumor_hvg/) { ds_name = "extended_tumor_hvg" }
        else if (input_file_path =~ /(?i)extended/) { ds_name = "extended" }
        
        log.info "Atlas input ${input_file_path} not found locally. Preparing automatic download & preprocessing for ${ds_name}..."
        PREPARE_DATASET(ds_name, file(params.download_script))
        atlas_input = PREPARE_DATASET.out
    } else {
        atlas_input = file(input_file_path)
    }
    
    TRAIN_SCANVI_ATLAS(file(params.train_script), PREPARE_ATLAS_CONFIG.out, atlas_input)
}

workflow SURGERY {
    setDefault('experiment_name', 'scanvi_surgery_smoke')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/scanvi_surgery/runs')
    setDefault('sif', '/data/containers/scvi-raytune-py313-cu12.sif')
    setDefault('bind_paths', '/datos,/data,/home,/datos/home')
    setDefault('train_script', "$baseDir/bin/run_surgery.py")
    setDefault('download_script', "$baseDir/bin/download_and_preprocess_dataset.py")
    setDefault('input_h5ad', '/datos/migccl/neto_maestria/luca_explore/surgeries/filtered_Trinks_Bishoff_2021_NSCLC.h5ad')
    setDefault('reference_model', '/datos/migccl/neto_maestria/luca_explore/LUCA_model/hvg_integrated_scvi_scanvi_tumor_model_b128_lay4_h1024_wuSep_epocs300-300')
    setDefault('slurm_cpus', 8)
    setDefault('slurm_memory', '64 GB')
    setDefault('slurm_time', '12:00:00')
    setDefault('slurm_gpus', 1)

    // ML defaults
    setDefault('seed', 0)
    setDefault('condition_key', 'dataset2')
    setDefault('cell_type_key', 'cell_type_tumor')
    setDefault('dataset_name', 'Bishoff_wu')
    setDefault('dataset_name_short', 'Bishoff')
    setDefault('dl_num_workers', 8)
    setDefault('float32_matmul_precision', 'high')
    setDefault('map_ensembl', true)
    setDefault('max_cells', null)
    setDefault('scanvi_train_params', [max_epochs: 200, batch_size: 128])
    setDefault('scanvi_early_stopping', true)
    setDefault('scanvi_early_stopping_kwargs', [early_stopping_monitor: 'elbo_validation', early_stopping_patience: 15, early_stopping_min_delta: 0.01])
    setDefault('scanvi_plan_kwargs', [reduce_lr_on_plateau: true, lr_patience: 9, lr_factor: 0.1])

    new File("${params.shared_run_root}/${params.experiment_name}").mkdirs()

    // Determine query dataset file (download if not exists and is a known dataset name)
    def query_file_path = params.input_h5ad
    def query_file
    def is_known_dataset = (query_file_path =~ /(?i)(zuani|deng|hu|bishoff|trinks)/)
    def file_exists = file(query_file_path).exists()
    
    if (!file_exists && is_known_dataset) {
        def ds_name = ""
        if (query_file_path =~ /(?i)zuani/) { ds_name = "Zuani" }
        else if (query_file_path =~ /(?i)deng/) { ds_name = "Deng" }
        else if (query_file_path =~ /(?i)(hu_zhang|hu2023|hu)/) { ds_name = "Hu" }
        else if (query_file_path =~ /(?i)(bishoff|trinks)/) { ds_name = "Bishoff" }
        
        log.info "Query dataset ${query_file_path} not found locally. Preparing automatic download & preprocessing for ${ds_name}..."
        PREPARE_DATASET(ds_name, file(params.download_script))
        query_file = PREPARE_DATASET.out
    } else {
        query_file = file(query_file_path)
    }

    PREPARE_SURGERY_CONFIG("scanvi_model")
    TRAIN_SCANVI_SURGERY(file(params.train_script), PREPARE_SURGERY_CONFIG.out, query_file, file(params.reference_model))
}

workflow RAYTUNE {
    setDefault('experiment_name', 'scvi_raytune_smoke')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/raytune_runs')
    setDefault('sif', '/data/containers/scvi-raytune-py313-cu12.sif')
    setDefault('bind_paths', '/datos,/data,/home,/datos/home')
    setDefault('train_script', "$baseDir/bin/train_scvi_raytune.py")
    setDefault('download_script', "$baseDir/bin/download_and_preprocess_dataset.py")
    setDefault('input_h5ad', '/data/luca_atlas/extended_tumor_hvg.h5ad')
    setDefault('slurm_cpus', 4)
    setDefault('slurm_memory', '32 GB')
    setDefault('slurm_time', '02:00:00')
    setDefault('slurm_gpus', 1)
    
    // defaults needed for raytune
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
    setDefault('search_space', [
        model_params: [
            n_hidden: [choice: [512, 1024, 2048]],
            n_layers: [choice: [5, 6, 7]],
            gene_likelihood: [choice: ['nb', 'zinb']]
        ],
        train_params: [
            batch_size: 128,
            plan_kwargs: [
                reduce_lr_on_plateau: [choice: [true, false]]
            ]
        ]
    ])
    setDefault('initial_points', [[
        model_params: [n_hidden: 2048, n_layers: 7, gene_likelihood: 'nb'],
        train_params: [batch_size: 128, plan_kwargs: [reduce_lr_on_plateau: true]]
    ]])
    setDefault('scheduler', 'asha')
    setDefault('searcher', 'hyperopt')
    setDefault('scheduler_kwargs', null)
    setDefault('searcher_kwargs', null)
    setDefault('save_checkpoints', true)
    setDefault('log_to_driver', true)
    setDefault('max_cells', null)
    setDefault('max_genes', null)
    setDefault('output_root', "${params.shared_run_root}/${params.experiment_name}/raytune_logs")

    new File("${params.shared_run_root}/${params.experiment_name}").mkdirs()

    // Determine raytune input dataset (download if not exists)
    def input_file_path = params.input_h5ad
    def raytune_input
    def is_atlas_dataset = (input_file_path =~ /(?i)(extended_tumor_hvg|extended)/)
    def file_exists = file(input_file_path).exists()
    
    if (!file_exists && is_atlas_dataset) {
        def ds_name = ""
        if (input_file_path =~ /(?i)extended_tumor_hvg/) { ds_name = "extended_tumor_hvg" }
        else if (input_file_path =~ /(?i)extended/) { ds_name = "extended" }
        
        log.info "Raytune input ${input_file_path} not found locally. Preparing automatic download & preprocessing for ${ds_name}..."
        PREPARE_DATASET(ds_name, file(params.download_script))
        raytune_input = PREPARE_DATASET.out
    } else {
        raytune_input = file(input_file_path)
    }

    RUN_SCVI_RAYTUNE(file(params.train_script), raytune_input)
}

workflow SUBCLUSTER {
    setDefault('experiment_name', 'luca_subcluster_default')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/subcluster/runs')
    subclusterDefaults()

    if (params.subcluster_input_h5ad && params.subcluster_cohort_manifest) {
        log.warn "Both subcluster_input_h5ad and subcluster_cohort_manifest are set; using subcluster_input_h5ad."
    }
    if (!params.subcluster_input_h5ad && !params.subcluster_cohort_manifest) {
        error "SUBCLUSTER requires subcluster_input_h5ad or subcluster_cohort_manifest"
    }
    if (params.subcluster_engine == 'scanpy' && params.subcluster_backend == 'gpu') {
        error "The scanpy subclustering engine only supports subcluster_backend=cpu"
    }

    def subclusterInput = file(params.subcluster_input_h5ad ?: params.subcluster_cohort_manifest)
    PREPARE_SUBCLUSTER_SCRIPT(file(params.subcluster_script))
    PREPARE_SUBCLUSTER_CONFIG(subclusterInput)
    RUN_SUBCLUSTER(
        PREPARE_SUBCLUSTER_SCRIPT.out.script_file,
        PREPARE_SUBCLUSTER_CONFIG.out.config_file,
        subclusterInput
    )
    APPLY_SUBCLUSTER_MAPPINGS(
        PREPARE_SUBCLUSTER_SCRIPT.out.script_file,
        RUN_SUBCLUSTER.out.raw_results,
        absoluteParamPath(params.early_mapping_json) ?: '',
        absoluteParamPath(params.late_mapping_json) ?: ''
    )
}

workflow GRAPH {
    setDefault('experiment_name', 'luca_graph_default')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/graph/runs')
    graphDefaults()
    if (!params.graph_input) {
        error "GRAPH requires graph_input (CSV or H5AD)"
    }
    def graphInput = file(params.graph_input)
    PREPARE_GRAPH_SCRIPT(file(params.graph_script))
    PREPARE_GRAPH_CONFIG(graphInput)
    RUN_GRAPH_PHASE(PREPARE_GRAPH_SCRIPT.out, PREPARE_GRAPH_CONFIG.out, graphInput)
}

workflow ANALYSIS {
    setDefault('experiment_name', 'luca_analysis_default')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/graph/runs')
    setDefault('graph_results_dir', null)
    setDefault('marker_results_dir', null)
    setDefault('circos_script', "$baseDir/bin/run_circos_plot.py")
    setDefault('cell_mappings_json', "$baseDir/metadata/cell_mappings.json")
    setDefault('run_circos', true)
    setDefault('cpdb_script', "$baseDir/bin/run_cellphonedb.py")
    setDefault('analysis_script', "$baseDir/bin/run_graph_analysis.py")
    setDefault('adata_dir', '/datos/migccl/neto_maestria/luca_explore/surgeries/')
    setDefault('cpdb_file', '/datos/migccl/neto_maestria/luca_explore/cellphonedb/cellphonedb_v4.0.0.zip')
    setDefault('marker_analysis_script', "$baseDir/bin/run_marker_analysis.py")
    setDefault('gsea_gmt', "$baseDir/metadata/h.all.v2023.2.Hs.symbols.gmt")
    setDefault('modal_de_script', "$baseDir/bin/modal_DE.py")
    setDefault('consensus_script', "$baseDir/bin/aggregate_consensus.py")

    if (!params.graph_results_dir) {
        error "ANALYSIS requires graph_results_dir"
    }
    if (!params.marker_results_dir) {
        error "ANALYSIS requires marker_results_dir"
    }

    def graphResults = file(params.graph_results_dir)
    def markerResults = file(params.marker_results_dir)

    def markersDir = file("${params.marker_results_dir}/final_auc_regions")

    PREPARE_MARKER_ANALYSIS_SCRIPT(file(params.marker_analysis_script))
    PREPARE_GSEA_GMT(file(params.gsea_gmt))
    PREPARE_GRAPH_ANALYSIS_SCRIPT(file(params.analysis_script))
    PREPARE_CELLPHONEDB_SCRIPT(file(params.cpdb_script))
    PREPARE_CIRCOS_SCRIPT(file(params.circos_script))
    PREPARE_CELL_MAPPINGS(file(params.cell_mappings_json))
    PREPARE_MODAL_DE_SCRIPT(file(params.modal_de_script))
    PREPARE_CONSENSUS_SCRIPT(file(params.consensus_script))

    RUN_CELLPHONEDB(
        PREPARE_CELLPHONEDB_SCRIPT.out, 
        markersDir, 
        markersDir, 
        graphResults, 
        params.adata_dir, 
        params.cpdb_file,
        PREPARE_MODAL_DE_SCRIPT.out
    )

    def aucFiles = Channel.fromPath("${params.marker_results_dir}/**/tumor-vs-all/*_auc.npy")
    RUN_MARKER_ANALYSIS(
        PREPARE_MARKER_ANALYSIS_SCRIPT.out,
        aucFiles,
        PREPARE_GSEA_GMT.out,
        PREPARE_MODAL_DE_SCRIPT.out
    )

    RUN_CONSENSUS_AGGREGATION(
        PREPARE_CONSENSUS_SCRIPT.out,
        markersDir,
        RUN_MARKER_ANALYSIS.out.results.collect()
    )

    RUN_GRAPH_ANALYSIS(PREPARE_GRAPH_ANALYSIS_SCRIPT.out, graphResults)

    if (params.run_circos) {
        RUN_CIRCOS_PLOT(PREPARE_CIRCOS_SCRIPT.out, graphResults, PREPARE_CELL_MAPPINGS.out)
    }

    // Collect all final outputs to the unified results/ directory structure
    def resultsDir = "${params.shared_run_root}/${params.experiment_name}/results"

    RUN_CELLPHONEDB.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/cellphonedb"
        cp -rn ${res}/. "${resultsDir}/cellphonedb/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    RUN_MARKER_ANALYSIS.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/marker_genes/analysis_plots"
        cp -rn ${res}/. "${resultsDir}/marker_genes/analysis_plots/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    RUN_GRAPH_ANALYSIS.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/graph_generation/analysis_plots"
        cp -rn ${res}/. "${resultsDir}/graph_generation/analysis_plots/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    if (params.run_circos) {
        RUN_CIRCOS_PLOT.out.results.map { res ->
            """
            mkdir -p "${resultsDir}/graph_generation/analysis_plots"
            cp -rn ${res}/. "${resultsDir}/graph_generation/analysis_plots/" || true
            """
        }.subscribe { cmd -> ['bash', '-c', cmd].execute() }
    }

    RUN_CONSENSUS_AGGREGATION.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/marker_genes/consensus"
        cp -r ${res}/. "${resultsDir}/marker_genes/consensus/"
        """
    }.subscribe { cmd -> 
        def proc = ['bash', '-c', cmd].execute()
        proc.waitFor()
    }
}

// ---------------------------------------------------------
// AUTO-CHAINED END-TO-END PIPELINE WORKFLOW
// ---------------------------------------------------------

workflow PIPELINE {
    // Shared execution defaults
    setDefault('experiment_name', 'scanvi_pipeline_run')
    setDefault('shared_run_root', '/datos/home/epaaso/slurm-gpu-jobs/raytune_runs')
    setDefault('sif', '/data/containers/scvi-raytune-py313-cu12.sif')
    setDefault('bind_paths', '/datos,/data,/home,/datos/home')
    
    // Scripts
    setDefault('raytune_script', "$baseDir/bin/train_scvi_raytune.py")
    setDefault('atlas_script', "$baseDir/bin/train_scanvi_atlas.py")
    setDefault('surgery_script', "$baseDir/bin/run_surgery.py")
    setDefault('apply_script', "$baseDir/bin/apply_raytune_best_config.py")
    setDefault('download_script', "$baseDir/bin/download_and_preprocess_dataset.py")
    setDefault('marker_script', "$baseDir/bin/run_marker_extraction.py")
    setDefault('cpdb_script', "$baseDir/bin/run_cellphonedb.py")
    setDefault('analysis_script', "$baseDir/bin/run_graph_analysis.py")
    setDefault('circos_script', "$baseDir/bin/run_circos_plot.py")
    setDefault('cell_mappings_json', "$baseDir/metadata/cell_mappings.json")
    setDefault('run_circos', true)
    setDefault('pipeline_cohort_manifest', "$baseDir/configs/subcluster_cohort_default.json")
    setDefault('adata_dir', '/datos/migccl/neto_maestria/luca_explore/surgeries/')
    setDefault('cpdb_file', '/datos/migccl/neto_maestria/luca_explore/cellphonedb/cellphonedb_v4.0.0.zip')
    setDefault('marker_analysis_script', "$baseDir/bin/run_marker_analysis.py")
    setDefault('gsea_gmt', "$baseDir/metadata/h.all.v2023.2.Hs.symbols.gmt")
    setDefault('modal_de_script', "$baseDir/bin/modal_DE.py")
    setDefault('consensus_script', "$baseDir/bin/aggregate_consensus.py")

    // Datasets
    setDefault('input_h5ad', '/data/luca_atlas/extended_tumor_hvg.h5ad')
    setDefault('query_h5ad', '/datos/migccl/neto_maestria/luca_explore/surgeries/filtered_Trinks_Bishoff_2021_NSCLC.h5ad')

    // The execution profile selects CPU/GPU backend and container.
    subclusterDefaults()
    graphDefaults()

    // Ray Tune parameters
    setDefault('metric', 'validation_loss')
    setDefault('mode', 'min')
    setDefault('num_samples', 50)
    setDefault('seed', 0)
    setDefault('batch_key', 'dataset')
    setDefault('labels_key', 'cell_type_tumor')
    setDefault('model', 'SCVI')
    setDefault('cpus_per_trial', 5)
    setDefault('gpus_per_trial', 0.25)
    setDefault('ray_cpus', 20)
    setDefault('ray_gpus', 1)
    setDefault('search_space', [
        model_params: [
            n_hidden: [choice: [512, 1024, 2048]],
            n_layers: [choice: [5, 6, 7]],
            gene_likelihood: [choice: ['nb', 'zinb']]
        ],
        train_params: [
            batch_size: 128,
            plan_kwargs: [
                reduce_lr_on_plateau: [choice: [true, false]]
            ]
        ]
    ])
    setDefault('initial_points', [[
        model_params: [n_hidden: 2048, n_layers: 7, gene_likelihood: 'nb'],
        train_params: [batch_size: 128, plan_kwargs: [reduce_lr_on_plateau: true]]
    ]])
    setDefault('scheduler', 'asha')
    setDefault('searcher', 'hyperopt')
    setDefault('scheduler_kwargs', null)
    setDefault('searcher_kwargs', null)
    setDefault('save_checkpoints', true)
    setDefault('log_to_driver', true)
    setDefault('max_cells', null)
    setDefault('max_genes', null)
    setDefault('output_root', "${params.shared_run_root}/${params.experiment_name}/raytune_logs")

    // Atlas parameters
    setDefault('unlabeled_category', 'Unknown')
    setDefault('origin_key', 'origin')
    setDefault('origin_values', ['tumor_primary'])
    setDefault('stage_key', 'uicc_stage')
    setDefault('stages', ['I', 'II', 'III', 'III or IV', 'IV'])
    setDefault('study_key', 'study')
    setDefault('exclude_studies', ['Goveia_Carmeliet_2020', 'Leader_Merad_2021', 'Guo_Zhang_2018', 'Wu_Zhou_2021'])
    setDefault('sample_key', 'sample')
    setDefault('split_batch_dataset', null)
    setDefault('counts_layer', 'count')
    setDefault('drop_layers', ['counts_length_scaled', 'count'])
    setDefault('drop_raw', true)
    setDefault('round_counts', false)
    setDefault('hvg_mode', 'compute')
    setDefault('hvg_n_top_genes', 6000)
    setDefault('hvg_flavor', 'seurat_v3')
    setDefault('hvg_batch_key', 'dataset')
    setDefault('hvg_column', 'is_highly_variable')
    setDefault('save_prepared_h5ad', false)
    setDefault('dl_num_workers', 8)
    setDefault('float32_matmul_precision', 'high')
    setDefault('numba_threads', 30)
    setDefault('compute_neighbors', true)
    setDefault('compute_leiden', true)
    setDefault('compute_umap', true)
    setDefault('n_neighbors', 8)
    setDefault('collapse_wu_dataset', true)
    setDefault('latent_obs_columns', ['dataset', 'study', 'sample', 'uicc_stage'])
    setDefault('tensorboard_logging', true)
    setDefault('tensorboard_log_dir', null)
    setDefault('reuse_scvi_model_if_exists', true)
    setDefault('scvi_trainer_kwargs', [gradient_clip_val: 10.0])
    setDefault('scanvi_trainer_kwargs', [gradient_clip_val: 10.0])

    // Surgery parameters
    setDefault('condition_key', 'dataset2')
    setDefault('cell_type_key', 'cell_type_tumor')
    setDefault('dataset_name', 'Bishoff_wu')
    setDefault('dataset_name_short', 'Bishoff')
    setDefault('map_ensembl', true)
    setDefault('scanvi_train_params', [max_epochs: 200, batch_size: 128])
    setDefault('scanvi_early_stopping', true)
    setDefault('scanvi_early_stopping_kwargs', [early_stopping_monitor: 'elbo_validation', early_stopping_patience: 15, early_stopping_min_delta: 0.01])
    setDefault('scanvi_plan_kwargs', [reduce_lr_on_plateau: true, lr_patience: 9, lr_factor: 0.1])

    // Execution resources
    setDefault('slurm_cpus', 8)
    setDefault('slurm_memory', '64 GB')
    setDefault('slurm_time', '24:00:00')
    setDefault('slurm_gpus', 1)

    new File("${params.shared_run_root}/${params.experiment_name}").mkdirs()

    // Determine raytune/atlas input dataset (download if not exists)
    def input_file_path = params.input_h5ad
    def atlas_input
    def is_atlas_dataset = (input_file_path =~ /(?i)(extended_tumor_hvg|extended)/)
    def file_exists = file(input_file_path).exists()
    
    if (!file_exists && is_atlas_dataset) {
        def ds_name = ""
        if (input_file_path =~ /(?i)extended_tumor_hvg/) { ds_name = "extended_tumor_hvg" }
        else if (input_file_path =~ /(?i)extended/) { ds_name = "extended" }
        
        log.info "Atlas input ${input_file_path} not found locally. Preparing automatic download & preprocessing for ${ds_name}..."
        PREPARE_DATASET(ds_name, file(params.download_script))
        atlas_input = PREPARE_DATASET.out
    } else {
        atlas_input = file(input_file_path)
    }

    // 1. Run Ray Tune hyperparameter exploration
    RUN_SCVI_RAYTUNE(file(params.raytune_script), atlas_input)
    
    // 2. Prepare the base Atlas configuration json
    PREPARE_ATLAS_CONFIG()
    
    // 3. Apply the best found config to the Atlas configuration json
    APPLY_RAYTUNE_CONFIG(PREPARE_ATLAS_CONFIG.out, RUN_SCVI_RAYTUNE.out.best_config_json, file(params.apply_script))
    
    // 4. Train the reference Atlas with the optimized hyperparameters
    TRAIN_SCANVI_ATLAS(file(params.atlas_script), APPLY_RAYTUNE_CONFIG.out, atlas_input)
    
    // 5. Fan out surgery over the configured cohort. If no cohort is supplied,
    // preserve query_h5ad as a single-query fallback.
    def cohortEntries = []
    if (params.pipeline_cohort_manifest) {
        def cohortSpec = new JsonSlurper().parse(file(params.pipeline_cohort_manifest).toFile())
        cohortEntries = cohortSpec.queries ?: []
    }
    if (!cohortEntries) {
        cohortEntries = [[
            name: params.dataset_name_short,
            input_h5ad: params.query_h5ad,
            download_name: params.dataset_name_short,
            batch_key: 'sample',
            stage_key: 'stage'
        ]]
    }
    def cohortChannel = Channel.fromList(cohortEntries)
    PREPARE_COHORT_QUERY(cohortChannel, file(params.download_script))
    RUN_COHORT_SURGERY(
        PREPARE_COHORT_QUERY.out,
        TRAIN_SCANVI_ATLAS.out.scanvi_model,
        file(params.surgery_script)
    )

    // 6. Build a runtime manifest from the newly generated atlas/query latent
    // representations and run early/late tumor-epithelial subclustering.
    PREPARE_PIPELINE_SUBCLUSTER_CONFIG(
        TRAIN_SCANVI_ATLAS.out.ref_latent_h5ad,
        RUN_COHORT_SURGERY.out.cohort_dir.collect()
    )
    PREPARE_SUBCLUSTER_SCRIPT(file(params.subcluster_script))
    RUN_SUBCLUSTER(
        PREPARE_SUBCLUSTER_SCRIPT.out.script_file,
        PREPARE_PIPELINE_SUBCLUSTER_CONFIG.out.config_file,
        PREPARE_PIPELINE_SUBCLUSTER_CONFIG.out.cohort_manifest
    )
    APPLY_SUBCLUSTER_MAPPINGS(
        PREPARE_SUBCLUSTER_SCRIPT.out.script_file,
        RUN_SUBCLUSTER.out.raw_results,
        absoluteParamPath(params.early_mapping_json) ?: '',
        absoluteParamPath(params.late_mapping_json) ?: ''
    )

    // 7. Build early/late global ARACNe networks, discover sample ecotypes,
    // and build ecotype-specific ARACNe networks.
    def adjustedCohort = APPLY_SUBCLUSTER_MAPPINGS.out.map { directory ->
        file("${directory}/query_latent_adjusted.h5ad")
    }
    PREPARE_GRAPH_SCRIPT(file(params.graph_script))
    PREPARE_GRAPH_CONFIG(adjustedCohort)
    RUN_GRAPH_PHASE(PREPARE_GRAPH_SCRIPT.out, PREPARE_GRAPH_CONFIG.out, adjustedCohort)
    
    // 8. Run Downstream Marker Extraction, Graph Analysis, and CellPhoneDB
    PREPARE_MARKER_EXTRACTION_SCRIPT(file(params.marker_script))
    PREPARE_MODAL_DE_SCRIPT(file(params.modal_de_script))

    marker_stages = Channel.of('early', 'late')
    marker_runs = RUN_MARKER_EXTRACTION(
        PREPARE_MARKER_EXTRACTION_SCRIPT.out, 
        RUN_GRAPH_PHASE.out.results, 
        params.adata_dir, 
        marker_stages,
        PREPARE_MODAL_DE_SCRIPT.out
    )
    early_markers = marker_runs.results.filter { it[0] == 'early' }.map { it[1] }
    late_markers = marker_runs.results.filter { it[0] == 'late' }.map { it[1] }

    // Gather all tumor-vs-all *_auc.npy files from RUN_MARKER_EXTRACTION
    def early_auc = early_markers.map { dir -> file("${dir}/tumor-vs-all/*_auc.npy") }.flatMap()
    def late_auc = late_markers.map { dir -> file("${dir}/tumor-vs-all/*_auc.npy") }.flatMap()
    def pipeline_auc = early_auc.concat(late_auc)

    PREPARE_MARKER_ANALYSIS_SCRIPT(file(params.marker_analysis_script))
    PREPARE_GSEA_GMT(file(params.gsea_gmt))
    PREPARE_GRAPH_ANALYSIS_SCRIPT(file(params.analysis_script))
    PREPARE_CELLPHONEDB_SCRIPT(file(params.cpdb_script))
    PREPARE_CIRCOS_SCRIPT(file(params.circos_script))
    PREPARE_CELL_MAPPINGS(file(params.cell_mappings_json))
    PREPARE_CONSENSUS_SCRIPT(file(params.consensus_script))

    RUN_MARKER_ANALYSIS(
        PREPARE_MARKER_ANALYSIS_SCRIPT.out,
        pipeline_auc,
        PREPARE_GSEA_GMT.out,
        PREPARE_MODAL_DE_SCRIPT.out
    )

    RUN_CONSENSUS_AGGREGATION(
        PREPARE_CONSENSUS_SCRIPT.out,
        early_markers.concat(late_markers).collect(),
        RUN_MARKER_ANALYSIS.out.results.collect()
    )

    RUN_CELLPHONEDB(
        PREPARE_CELLPHONEDB_SCRIPT.out, 
        early_markers, 
        late_markers, 
        RUN_GRAPH_PHASE.out.results, 
        params.adata_dir, 
        params.cpdb_file,
        PREPARE_MODAL_DE_SCRIPT.out
    )
    RUN_GRAPH_ANALYSIS(PREPARE_GRAPH_ANALYSIS_SCRIPT.out, RUN_GRAPH_PHASE.out.results)

    if (params.run_circos) {
        RUN_CIRCOS_PLOT(PREPARE_CIRCOS_SCRIPT.out, RUN_GRAPH_PHASE.out.results, PREPARE_CELL_MAPPINGS.out)
    }

    // Collect all final outputs to the unified results/ directory structure
    def resultsDir = "${params.shared_run_root}/${params.experiment_name}/results"
    
    RUN_GRAPH_PHASE.out.results.map { graph_results ->
        """
        set -euo pipefail
        mkdir -p "${resultsDir}/graph_generation" "${resultsDir}/ecotype_clustering"
        cp -rn ${graph_results}/global/. "${resultsDir}/graph_generation/" || true
        cp -rn ${graph_results}/ecotype_graphs/. "${resultsDir}/graph_generation/" || true
        cp -rn ${graph_results}/ecotype/. "${resultsDir}/ecotype_clustering/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }
    
    early_markers.map { res ->
        """
        mkdir -p "${resultsDir}/marker_genes"
        cp -rn ${res}/. "${resultsDir}/marker_genes/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    late_markers.map { res ->
        """
        mkdir -p "${resultsDir}/marker_genes"
        cp -rn ${res}/. "${resultsDir}/marker_genes/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    RUN_CELLPHONEDB.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/cellphonedb"
        cp -rn ${res}/. "${resultsDir}/cellphonedb/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    RUN_MARKER_ANALYSIS.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/marker_genes/analysis_plots"
        cp -rn ${res}/. "${resultsDir}/marker_genes/analysis_plots/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    RUN_GRAPH_ANALYSIS.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/graph_generation/analysis_plots"
        cp -rn ${res}/. "${resultsDir}/graph_generation/analysis_plots/" || true
        """
    }.subscribe { cmd -> ['bash', '-c', cmd].execute() }

    if (params.run_circos) {
        RUN_CIRCOS_PLOT.out.results.map { res ->
            """
            mkdir -p "${resultsDir}/graph_generation/analysis_plots"
            cp -rn ${res}/. "${resultsDir}/graph_generation/analysis_plots/" || true
            """
        }.subscribe { cmd -> ['bash', '-c', cmd].execute() }
    }

    RUN_CONSENSUS_AGGREGATION.out.results.map { res ->
        """
        mkdir -p "${resultsDir}/marker_genes/consensus"
        cp -r ${res}/. "${resultsDir}/marker_genes/consensus/"
        """
    }.subscribe { cmd -> 
        def proc = ['bash', '-c', cmd].execute()
        proc.waitFor()
    }
}

// Default entry point runs the chained pipeline
workflow {
    PIPELINE()
}
