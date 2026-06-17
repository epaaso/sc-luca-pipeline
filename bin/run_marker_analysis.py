#!/usr/bin/env python3
"""Run marker genes plotting and GSEA analysis from a pre-calculated .npy marker file."""

import argparse
import logging
import os
import sys
import numpy as np

# Ensure the script can import modal_DE
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modal_DE import CommonConfig, VisualizerConfig, DEVisualizer

def main() -> int:
    parser = argparse.ArgumentParser(description="Plot marker genes and run GSEA on .npy markers.")
    parser.add_argument("--marker-file", required=True, help="Path to .npy marker file")
    parser.add_argument("--gsea-gmt", required=True, help="Path to GSEA GMT file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--name", required=True, help="Name/Prefix for the output files")
    parser.add_argument("--time", required=True, help="Time/Stage label (e.g. I-II or III-IV)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    logging.info("Starting marker analysis for %s using %s", args.marker_file, args.gsea_gmt)

    if not os.path.exists(args.marker_file):
        logging.error("Marker file not found: %s", args.marker_file)
        return 1

    # Load pre-calculated marker data
    try:
        de_region = np.load(args.marker_file, allow_pickle=True).item()
    except Exception as exc:
        logging.error("Failed to load numpy file %s: %s", args.marker_file, exc)
        return 1

    # Extract available cell types from structured array
    try:
        types = list(de_region["scores"].dtype.names)
    except KeyError:
        logging.error("The .npy file does not have the expected structured 'scores' key.")
        return 1

    if not types:
        logging.warning("No cell types found in %s, skipping analysis.", args.marker_file)
        return 0

    os.makedirs(args.output_dir, exist_ok=True)

    # Initialize configuration
    common = CommonConfig(
        ext_name=args.name,
        name=args.name,
        time=args.time,
        w_folder=args.output_dir
    )

    vis_config = VisualizerConfig(
        region_mapping="",
        load_gsea=False,
        load_gsea_heatmap=False,
        gsea_gmt=args.gsea_gmt
    )

    visualizer = DEVisualizer(common, vis_config)

    # Run GSEA and marker gene plotting
    logging.info("Plotting marker genes for cell types: %s", types)
    try:
        visualizer.plot_marker_genes(de_region, valid_types=types)
    except Exception as exc:
        logging.error("Error plotting marker genes: %s", exc)

    logging.info("Running and plotting GSEA...")
    try:
        visualizer.plot_gsea(de_region, valid_types=types)
    except Exception as exc:
        logging.error("Error plotting GSEA: %s", exc)

    logging.info("Marker analysis completed successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
