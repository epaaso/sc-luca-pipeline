#!/usr/bin/env python3
"""Run circos blockmodel plots using graph-tool."""

import argparse
import json
import logging
import os
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# This script is strictly for a graph-tool environment
import graph_tool.all as gt

LOGGER = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph-file", type=Path, required=True, help="Path to Pearson-annotated network CSV")
    parser.add_argument("--metadata", type=Path, help="Path to cell_mappings.json for colors")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--time-label", type=str, required=True, help="Time/Stage label for titles")
    parser.add_argument("--mi-threshold", type=float, default=0.0, help="MI threshold for plotting")
    return parser.parse_args()

def cell_category_mapping(cell_type: str, metadata: dict) -> str:
    if "Tumor" in cell_type:
        return "tumoral"
    cat_map = metadata.get("cell_categories", {})
    mapping = {ct: cat for cat, cts in cat_map.items() for ct in cts}
    return mapping.get(cell_type, "unknown")

def process_blockmodel(graph_file: Path, output_dir: Path, time_label: str, mi_threshold: float, metadata: dict):
    LOGGER.info("Generating blockmodel for %s", graph_file)
    try:
        g = gt.load_graph_from_csv(
            str(graph_file),
            csv_options={"delimiter": "\t"},
            eprop_types=["float", "float", "float", "float", "string"],
            eprop_names=["MI", "pvals", "pearson", "pvalue", "sign"],
        )
    except Exception as e:
        LOGGER.error("Failed to load graph with graph-tool: %s", e)
        return

    edges = list(g.edges())
    seen = set()
    for e in edges:
        src, dst = int(e.source()), int(e.target())
        if (src, dst) in seen or (dst, src) in seen:
            g.remove_edge(e)
        else:
            seen.add((src, dst))

    state = gt.minimize_nested_blockmodel_dl(g)
    gt.mcmc_anneal(state, beta_range=(1, 20), niter=1500, mcmc_equilibrate_args=dict(force_niter=20))

    color_map = metadata.get("color_map", {})
    color_map.setdefault("unknown", "#bcbcbc")

    g.vp["cat"] = g.new_vertex_property("string")
    g.vp["color"] = g.new_vertex_property("string")
    g.vp["label"] = g.new_vertex_property("string")
    for v in g.vertices():
        name = g.vp.name[v] if "name" in g.vp else str(int(v))
        category = cell_category_mapping(name, metadata)
        color = color_map.get(category, color_map["unknown"])
        g.vp.cat[v] = category
        g.vp.color[v] = color
        g.vp.label[v] = name
    
    edge_width = g.new_edge_property("float")
    width_divisor = 8.0
    if "MI" in g.ep:
        adjusted = np.maximum(g.ep.MI.a - mi_threshold, 0.0)
        edge_width.a = adjusted / width_divisor
    else:
        edge_width.a = np.full(g.num_edges(), 1.0 / width_divisor)

    edge_dash = g.new_edge_property("vector<float>")
    if "sign" in g.ep:
        for e in g.edges():
            sign = g.ep.sign[e]
            if sign == "positive_significant":
                edge_dash[e] = [3, 3, 3]
            elif sign == "negative_significant":
                edge_dash[e] = [10, 10, 10]
            else:
                edge_dash[e] = [1, 0, 0]

    edge_gradient = g.new_edge_property("vector<float>")
    for e in g.edges():
        src = e.source()
        rgba = list(mcolors.to_rgba(g.vp.color[src]))
        edge_gradient[e] = [0, rgba[0], rgba[1], rgba[2], 0.7, 1, rgba[0], rgba[1], rgba[2], 0.2]

    fig, ax = plt.subplots(figsize=(12, 10))
    drawing = state.draw(
        vertex_fill_color=g.vp.color,
        edge_pen_width=edge_width,
        vertex_text=g.vp.label,
        vertex_font_size=0.05,
        edge_dash_style=edge_dash,
        vertex_text_position="centered",
        edge_gradient=edge_gradient,
        vertex_pen_width=0.0,
        mplfig=ax,
    )
    if drawing and len(drawing) > 0:
        drawing[0].fit_view(margin=0.25, yflip=True)
    ax.set_title(time_label)
    ax.set_axis_off()

    handles = []
    for category, color in color_map.items():
        if category == "unknown":
            continue
        handles.append(mpatches.Patch(color=color, label=category))
    if "sign" in g.ep:
        handles.append(Line2D([], [], color="black", linestyle="--", label="negative_significant"))
    handles.append(Line2D([], [], color="black", linestyle="dotted", label="positive_significant"))
    if handles:
        ax.legend(handles=handles, bbox_to_anchor=(1, 0.5), loc="upper left")
    plt.subplots_adjust(right=0.8, left=0.01, bottom=0.01, top=0.95)

    out_file = output_dir / f"circos_blockmodel_{time_label}.png"
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    LOGGER.info("Saved blockmodel to %s", out_file)

def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    metadata = {}
    if args.metadata and args.metadata.exists():
        with args.metadata.open() as fh:
            metadata = json.load(fh)

    process_blockmodel(args.graph_file, args.output_dir, args.time_label, args.mi_threshold, metadata)

if __name__ == "__main__":
    main()
