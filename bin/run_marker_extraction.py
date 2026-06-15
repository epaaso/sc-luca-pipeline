#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import json
import pathlib
import sys
from collections.abc import Mapping
from typing import Any
import yaml

from modal_DE import MARKER_CONTRASTS, _get_de_impl, common_kwargs

def read_structured(path: str | pathlib.Path) -> dict[str, Any]:
    path = pathlib.Path(path)
    with path.open() as handle:
        return yaml.safe_load(handle) if path.suffix.lower() in {".yaml", ".yml"} else json.load(handle)

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Pipeline config file")
    parser.add_argument("--adata-dir", required=True, type=pathlib.Path, help="Directory containing filtered H5AD files")
    parser.add_argument("--membership", required=True, type=pathlib.Path, help="Ecotype membership CSV")
    parser.add_argument("--stage", required=True, choices=["early", "late"], help="Stage (early/late)")
    parser.add_argument("--output_dir", required=True, type=pathlib.Path, help="Output directory")
    args = parser.parse_args()

    cfg = read_structured(args.config)
    output_root = args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)

    time_mapping = {"early": "I-II", "late": "III-IV"}
    time = time_mapping[args.stage]

    # Read membership to find all (dataset, cluster)
    import pandas as pd
    df = pd.read_csv(args.membership)
    # dataset column might not be exactly named, use the one from config
    dataset_col = cfg.get("dataset_key", "dataset")
    if dataset_col not in df.columns:
        print(f"Warning: {dataset_col} not found in membership. Using first dataset column.")
        dataset_col = [c for c in df.columns if c.lower() in ("dataset", "study", "batch")][0]
    
    unique_runs = df[[dataset_col, "membership"]].drop_duplicates()
    
    contrasts = ["tumor-vs-all", "normal-vs-normal"]
    
    for _, row in unique_runs.iterrows():
        dataset = str(row[dataset_col])
        cluster = int(row["membership"])
        
        # Determine file_ext_name
        # Typically the dataset name in the file might be prefixed.
        # We will assume it matches exactly for now.
        
        for contrast in contrasts:
            kwargs = {**common_kwargs}
            kwargs.update({
                "ext_name": dataset,
                "name": dataset,
                "time": time,
                "cell_key": cfg.get("cell_key", "cell_type_adjusted"),
                "stage_key": cfg.get("stage_key", "stage"),
                "pred_name": dataset,
                "obs_unique": False,
                "obs_has_name": False,
                "contrast": contrast,
                "w_folder": str(output_root),
                "cluster_id": cluster,
                "membership_csv": str(args.membership),
                "log_layer": "do_log1p",
                "load_pair": False,
                "skip_visualization": True,
                "backup_dir": str(args.adata_dir),
                "file_ext_name": dataset
            })

            if contrast == "normal-vs-normal":
                kwargs.update(parallel_pair=True, n_jobs_inner=1, num_processes=4, max_cells_per_type=2000)

            try:
                _get_de_impl(**kwargs)
                gc.collect()
            except Exception as error:
                print(f"ERROR running DE for {dataset} cluster {cluster} {contrast}: {error}")
                manifest = output_root / contrast / f"{time}_{dataset}_cluster_{cluster}_manifest.json"
                manifest.parent.mkdir(parents=True, exist_ok=True)
                manifest.write_text(json.dumps({
                    "dataset": dataset,
                    "stage": time,
                    "cluster": cluster,
                    "contrast": contrast,
                    "status": "failed",
                    "error": str(error),
                }, indent=2))

    return 0

if __name__ == "__main__":
    sys.exit(main())
