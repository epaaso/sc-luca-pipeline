#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import pathlib
import random
import shutil
import sys
import time
from collections.abc import Mapping, Sequence
from importlib.metadata import version
from typing import Any

import pandas as pd
import ray
import scanpy as sc
import scvi
import torch
import yaml
from ray import tune
from scvi import autotune


def read_config(path: pathlib.Path) -> dict[str, Any]:
    with path.open() as handle:
        if path.suffix.lower() in {".yaml", ".yml"}:
            return yaml.safe_load(handle)
        return json.load(handle)


def as_jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, pathlib.Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(k): as_jsonable(v) for k, v in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [as_jsonable(v) for v in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return str(value)


def flatten_dict(data: Mapping[str, Any], prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in data.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, Mapping):
            out.update(flatten_dict(value, name))
        else:
            out[name] = as_jsonable(value)
    return out


def tune_space(value: Any) -> Any:
    if isinstance(value, Mapping):
        keys = set(value.keys())
        if keys == {"choice"}:
            return tune.choice(value["choice"])
        if keys == {"grid_search"}:
            return tune.grid_search(value["grid_search"])
        if keys == {"uniform"}:
            low, high = value["uniform"]
            return tune.uniform(low, high)
        if keys == {"loguniform"}:
            low, high = value["loguniform"]
            return tune.loguniform(low, high)
        if keys == {"randint"}:
            low, high = value["randint"]
            return tune.randint(low, high)
        return {k: tune_space(v) for k, v in value.items()}
    if isinstance(value, list):
        return [tune_space(v) for v in value]
    return value


def write_json(path: pathlib.Path, data: Any) -> None:
    path.write_text(json.dumps(as_jsonable(data), indent=2, sort_keys=True))


def result_rows(result_grid: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in result_grid:
        row: dict[str, Any] = {}
        row["trial_path"] = getattr(result, "path", None)
        row.update({f"metric.{k}": as_jsonable(v) for k, v in result.metrics.items()})
        row.update({f"config.{k}": v for k, v in flatten_dict(result.config).items()})
        rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=pathlib.Path)
    args = parser.parse_args()

    cfg = read_config(args.config)
    cwd = pathlib.Path.cwd()
    experiment_name = cfg["experiment_name"]
    run_dir = pathlib.Path(cfg["run_dir"])
    run_dir.mkdir(parents=True, exist_ok=True)
    local_ray_tmp = pathlib.Path(os.environ.get("TMPDIR", "/tmp")) / f"ray-{experiment_name}-{os.getpid()}"
    local_ray_tmp.mkdir(parents=True, exist_ok=True)

    seed = int(cfg.get("seed", 0))
    random.seed(seed)
    scvi.settings.seed = seed
    torch.set_float32_matmul_precision("high")

    versions = {
        "python": sys.version,
        "anndata": version("anndata"),
        "pandas": pd.__version__,
        "ray": ray.__version__,
        "scanpy": version("scanpy"),
        "scib_metrics": version("scib-metrics"),
        "scvi": scvi.__version__,
        "torch": torch.__version__,
        "torch_cuda_available": torch.cuda.is_available(),
        "torch_cuda_device_count": torch.cuda.device_count(),
        "torch_cuda_device_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
    }
    write_json(cwd / "versions.json", versions)
    write_json(run_dir / "versions.json", versions)

    adata = sc.read_h5ad(cfg["input_h5ad"])
    max_cells = cfg.get("max_cells")
    if max_cells and adata.n_obs > int(max_cells):
        sc.pp.subsample(adata, n_obs=int(max_cells), random_state=seed, copy=False)
    max_genes = cfg.get("max_genes")
    if max_genes and adata.n_vars > int(max_genes):
        adata = adata[:, : int(max_genes)].copy()

    model_name = cfg.get("model", "SCVI")
    if model_name != "SCVI":
        raise ValueError(f"Only model='SCVI' is implemented, got {model_name!r}")
    model_cls = scvi.model.SCVI

    setup_kwargs: dict[str, Any] = {"batch_key": cfg["batch_key"]}
    labels_key = cfg.get("labels_key")
    if labels_key:
        setup_kwargs["labels_key"] = labels_key
    model_cls.setup_anndata(adata, **setup_kwargs)

    search_space = tune_space(cfg["search_space"])
    searcher_kwargs = dict(cfg.get("searcher_kwargs") or {})
    initial_points = cfg.get("initial_points") or []
    if initial_points:
        searcher_kwargs["points_to_evaluate"] = initial_points

    ray_cfg = cfg.get("ray") or {}
    ray.init(
        ignore_reinit_error=True,
        include_dashboard=False,
        log_to_driver=bool(cfg.get("log_to_driver", False)),
        num_cpus=int(ray_cfg.get("num_cpus") or os.cpu_count() or 1),
        num_gpus=float(ray_cfg.get("num_gpus", 0)),
        _temp_dir=str(local_ray_tmp),
    )

    started = time.time()
    try:
        scvi_tuner = autotune.run_autotune(
            model_cls,
            adata,
            metrics=cfg.get("metric", "validation_loss"),
            mode=cfg.get("mode", "min"),
            search_space=search_space,
            num_samples=int(cfg.get("num_samples", 1)),
            scheduler=cfg.get("scheduler", "asha"),
            searcher=cfg.get("searcher", "hyperopt"),
            seed=seed,
            resources=cfg.get("resources") or {"cpu": 1, "gpu": 0},
            experiment_name=experiment_name,
            logging_dir=str(run_dir),
            save_checkpoints=bool(cfg.get("save_checkpoints", False)),
            scheduler_kwargs=cfg.get("scheduler_kwargs") or None,
            searcher_kwargs=searcher_kwargs or None,
            log_to_driver=bool(cfg.get("log_to_driver", False)),
            ignore_reinit_error=True,
        )
    finally:
        ray.shutdown()

    result_grid = scvi_tuner.result_grid
    rows = result_rows(result_grid)
    summary = pd.DataFrame(rows)
    summary.insert(0, "experiment_name", experiment_name)
    summary.insert(1, "elapsed_seconds_total", time.time() - started)
    summary.to_csv(cwd / "summary.csv", index=False)
    summary.to_csv(run_dir / "summary.csv", index=False)

    metric = cfg.get("metric", "validation_loss")
    mode = cfg.get("mode", "min")
    best = result_grid.get_best_result(metric=metric, mode=mode)
    write_json(cwd / "best_config.json", best.config)
    write_json(run_dir / "best_config.json", best.config)
    write_json(cwd / "best_result.json", {"metrics": best.metrics, "path": getattr(best, "path", None)})
    write_json(run_dir / "best_result.json", {"metrics": best.metrics, "path": getattr(best, "path", None)})

    ray_log_dir_file = cwd / "ray_log_dir.txt"
    ray_log_dir_file.write_text(str(run_dir) + "\n")
    shutil.copy2(ray_log_dir_file, run_dir / "ray_log_dir.txt")
    print(f"Finished {experiment_name}; results in {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
