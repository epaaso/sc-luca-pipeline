#!/usr/bin/env python3
"""Run graph analysis to produce normal (linear) and semi-log degree distribution plots."""

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx

LOGGER = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph-file", type=Path, required=True, help="Path to Pearson-annotated network CSV")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--time-label", type=str, required=True, help="Time/Stage label for titles")
    return parser.parse_args()

def plot_degree_distribution(graph_file: Path, output_dir: Path, title: str):
    LOGGER.info("Plotting degree distribution for %s", graph_file)
    try:
        df = pd.read_csv(graph_file, sep="\t")
        if df.empty:
            LOGGER.warning("Graph file %s is empty, skipping plot", graph_file)
            return

        # Handle header or no-header input
        if "source" in df.columns:
            sources = df["source"].astype(str)
            targets = df["target"].astype(str)
        else:
            sources = df.iloc[:, 0].astype(str)
            targets = df.iloc[:, 1].astype(str)
            
        G = nx.Graph()
        for s, t in zip(sources, targets):
            G.add_edge(s, t)
            
        degrees = [d for n, d in G.degree()]
        if not degrees:
            return
        degree_counts = np.bincount(degrees)
        degrees_unique = np.nonzero(degree_counts)[0]
        counts = degree_counts[degrees_unique]

        # 1. Normal (linear-linear) plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(degrees_unique, counts, alpha=0.7)
        ax.set_xlabel('Degree (linear scale)')
        ax.set_ylabel('Count (linear scale)')
        ax.set_title(f'Degree Distribution (Normal) - {title}')
        ax.grid(True, which="both", ls="--", alpha=0.2)
        normal_file = output_dir / f"degree_distribution_normal_{title}.png"
        fig.savefig(normal_file, dpi=200, bbox_inches='tight')
        plt.close(fig)
        LOGGER.info("Saved normal plot to %s", normal_file)

        # 2. Semi-log plot (x linear, y log)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(degrees_unique, counts, alpha=0.7)
        ax.set_yscale('log')
        ax.set_xlabel('Degree (linear scale)')
        ax.set_ylabel('Count (log scale)')
        ax.set_title(f'Degree Distribution (Semi-log) - {title}')
        ax.grid(True, which="both", ls="--", alpha=0.2)
        semilog_file = output_dir / f"degree_distribution_semilog_{title}.png"
        fig.savefig(semilog_file, dpi=200, bbox_inches='tight')
        plt.close(fig)
        LOGGER.info("Saved semi-log plot to %s", semilog_file)

    except Exception as exc:
        LOGGER.warning("Could not plot degree distribution: %s", exc)

def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_degree_distribution(args.graph_file, args.output_dir, args.time_label)

if __name__ == "__main__":
    main()
