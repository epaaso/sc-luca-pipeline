# Pipeline Skills

Use this file as procedural guidance for agents working on the unified Nextflow pipeline.

## Skill: Monitor Active Runs

Goal: Determine whether the current run is healthy, stuck, failed, or done.

The execution can be local or on a remote GPU. To inspect the actual python process, you must check the host where the compute is taking place.

### For Remote Slurm Jobs
```bash
ssh epaaso@inmegen-gpu '
local=/data/scanvi_atlas_runs/<experiment>
date
ps -p "$(cat "$local/direct_train.pid")" -o pid,stat,etime,cmd || true
pgrep -af "train_scanvi_atlas.py|run_surgery.py|train_scvi_raytune.py|tensorboard.main" || true
nvidia-smi --query-gpu=timestamp,name,utilization.gpu,memory.used,power.draw --format=csv,noheader
tail -120 "$local/train.log" 2>/dev/null || true
ls -lh "$local"/metrics.json "$local"/*.h5ad "$local"/scanvi_model/model.pt 2>/dev/null || true
'
```

### For Local Jobs
```bash
# Check running python processes related to the pipeline locally
pgrep -af "train_scanvi_atlas.py|run_surgery.py|train_scvi_raytune.py|tensorboard.main" || true
# If you have a local GPU:
nvidia-smi --query-gpu=timestamp,name,utilization.gpu,memory.used,power.draw --format=csv,noheader
# Monitor logs in the run directory configured via `shared_run_root`
```

Healthy signs:
- The script is alive and consuming CPU/GPU.
- `train.log` shows increasing SCVI or SCANVI epochs.
- Loss is finite and not exploding.

Failure signs:
- `train.log` ends in a Python traceback.
- Process is dead and final outputs are missing.

## Skill: Check TensorBoard

Goal: Verify native TensorBoard logging is available and exposed.

TensorBoard logs are usually output to `tensorboard_native/` inside the run directory.

### Remote Access
```bash
ssh -f -N -L 127.0.0.1:6007:127.0.0.1:6007 epaaso@inmegen-gpu
curl -sS --max-time 3 -I http://127.0.0.1:6007/ | head -1
```

### Local Access
If the run is local, you can launch TensorBoard directly on the run directory:
```bash
tensorboard --logdir <run_dir>/tensorboard_native --port 6007
```

## Skill: Fix Too Many Open Files

Goal: Resolve dataloader worker failures during PyTorch multiprocessing.

If a run fails with `RuntimeError: Too many open files. Communication with the workers is no longer possible.`:
1. Limit workers by configuring `dl_num_workers: 2` in the config.
2. Set `torch_multiprocessing_sharing_strategy: file_system` in the config.
3. Ensure `ulimit -n 65535` is set in the execution environment if possible.

## Skill: Fix Numeric Config Values

Goal: Prevent PyTorch/scvi optimizer errors caused by YAML/JSON string numbers.

If you encounter `TypeError: '<=' not supported between instances of 'float' and 'str'`, ensure the python scripts parse string numbers. The scripts (`run_surgery.py`, `train_scanvi_atlas.py`) include `coerce_numeric_strings()` to cast `"1e-4"` into a float before passing kwargs to scvi. Check the `run_config.json` payload in the run directory.
