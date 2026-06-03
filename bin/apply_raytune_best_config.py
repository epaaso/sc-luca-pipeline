#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import pathlib
from typing import Any

import yaml


def read_yaml(path: pathlib.Path) -> dict[str, Any]:
    with path.open() as handle:
        data = yaml.safe_load(handle)
    return data or {}


def write_yaml(path: pathlib.Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False))


def read_json(path: pathlib.Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def merge_dict(target: dict[str, Any], source: dict[str, Any]) -> dict[str, Any]:
    out = dict(target)
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge_dict(out[key], value)
        else:
            out[key] = value
    return out


def apply_best_config(
    base_config: dict[str, Any],
    best_config: dict[str, Any],
    *,
    copy_max_epochs: bool,
) -> dict[str, Any]:
    out = copy.deepcopy(base_config)

    model_params = best_config.get("model_params") or {}
    out["scvi_model_params"] = merge_dict(out.get("scvi_model_params") or {}, model_params)

    train_params = dict(best_config.get("train_params") or {})
    plan_kwargs = train_params.pop("plan_kwargs", None)
    if not copy_max_epochs:
        train_params.pop("max_epochs", None)
    out["scvi_train_params"] = merge_dict(out.get("scvi_train_params") or {}, train_params)
    if plan_kwargs:
        out["scvi_plan_kwargs"] = merge_dict(out.get("scvi_plan_kwargs") or {}, plan_kwargs)

    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a scanvi_atlas params YAML from a scvi_raytune best_config.json."
    )
    parser.add_argument("--base", required=True, type=pathlib.Path, help="Base scanvi_atlas YAML config.")
    parser.add_argument("--best-config", required=True, type=pathlib.Path, help="Ray Tune best_config.json.")
    parser.add_argument("--output", required=True, type=pathlib.Path, help="Output scanvi_atlas YAML config.")
    parser.add_argument("--experiment-name", help="Optional experiment_name for the output config.")
    parser.add_argument(
        "--copy-max-epochs",
        action="store_true",
        help="Also copy train_params.max_epochs from Ray Tune into scvi_train_params.",
    )
    args = parser.parse_args()

    config = apply_best_config(
        read_yaml(args.base),
        read_json(args.best_config),
        copy_max_epochs=args.copy_max_epochs,
    )
    if args.experiment_name:
        config["experiment_name"] = args.experiment_name
    write_yaml(args.output, config)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
