#!/usr/bin/env python3
"""Run graph analysis to produce degree distribution plots, cell type coabundance plots, and sample PCA plots."""

import argparse
import logging
import re
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
    parser.add_argument("--groups-file", type=Path, default=None, help="Path to groups_*.csv file")
    parser.add_argument("--membership-file", type=Path, default=None, help="Path to membership_*.csv file")
    parser.add_argument("--abundance-file", type=Path, default=None, help="Path to abundance_*.csv file")
    parser.add_argument("--global-membership-file", type=Path, default=None, help="Path to membership_*.csv file for global PCA")
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

def plot_pair_helper(df: pd.DataFrame, cell1: str, cell2: str, output_path: Path, dataset_map: pd.Series, axis_label: str):
    available_cells = set(df.columns)
    if cell1 not in available_cells or cell2 not in available_cells:
        LOGGER.warning("Cell '%s' or '%s' not in available columns", cell1, cell2)
        return

    x = df[cell1]
    y = df[cell2]
    
    # Prepare colors
    datasets = dataset_map.reindex(df.index).fillna("unknown")
    unique_datasets = sorted(datasets.unique())
    cmap = plt.get_cmap("tab10", len(unique_datasets) if len(unique_datasets) > 0 else 1)
    color_map = {d: cmap(i) for i, d in enumerate(unique_datasets)}
    colors = datasets.map(color_map)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    ((ax1, ax2), (ax3, ax4)) = axes
    
    # Plot 1: Log Y, Linear X
    ax1.scatter(x, y, c=colors, alpha=0.7, s=30, edgecolor='k', linewidth=0.3)
    ax1.set_yscale('log')
    ax1.set_xscale('linear')
    ax1.set_xlabel(f"{cell1} ({axis_label})")
    ax1.set_ylabel(f"{cell2} ({axis_label}) [Log]")
    ax1.set_title(f"{cell1} vs {cell2}\n(Log Y-axis)")
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    
    # Plot 2: Log X, Linear Y
    ax2.scatter(x, y, c=colors, alpha=0.7, s=30, edgecolor='k', linewidth=0.3)
    ax2.set_xscale('log')
    ax2.set_yscale('linear')
    ax2.set_xlabel(f"{cell1} ({axis_label}) [Log]")
    ax2.set_ylabel(f"{cell2} ({axis_label})")
    ax2.set_title(f"{cell1} vs {cell2}\n(Log X-axis)")
    ax2.grid(True, which="both", ls="-", alpha=0.2)

    # Plot 3: Linear X, Linear Y
    ax3.scatter(x, y, c=colors, alpha=0.7, s=30, edgecolor='k', linewidth=0.3)
    ax3.set_xscale('linear')
    ax3.set_yscale('linear')
    ax3.set_xlabel(f"{cell1} ({axis_label})")
    ax3.set_ylabel(f"{cell2} ({axis_label})")
    ax3.set_title(f"{cell1} vs {cell2}\n(Linear)")
    ax3.grid(True, which="both", ls="-", alpha=0.2)

    # Plot 4: Log X, Log Y
    ax4.scatter(x, y, c=colors, alpha=0.7, s=30, edgecolor='k', linewidth=0.3)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel(f"{cell1} ({axis_label}) [Log]")
    ax4.set_ylabel(f"{cell2} ({axis_label}) [Log]")
    ax4.set_title(f"{cell1} vs {cell2}\n(Log-Log)")
    ax4.grid(True, which="both", ls="-", alpha=0.2)
    
    # Add legend
    if unique_datasets:
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[d], label=d, markersize=10) for d in unique_datasets]
        fig.legend(handles=handles, loc='center left', title="Dataset", bbox_to_anchor=(0.85, 0.5))
        plt.tight_layout(rect=[0, 0, 0.85, 1])
    else:
        plt.tight_layout()
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    plt.close(fig)
    LOGGER.info("Saved coabundance plot to %s", output_path)

def plot_cell_type_coabundance(graph_file: Path, groups_file: Path, membership_file: Path, output_dir: Path, time_label: str):
    cluster_match = re.search(r'cluster_?(\d+)', time_label)
    if not cluster_match:
        LOGGER.info("No cluster ID found in time-label %s, skipping coabundance plots", time_label)
        return
    
    cluster_id = int(cluster_match.group(1))
    
    LOGGER.info("Running coabundance plotting for cluster=%d in %s", cluster_id, graph_file)
    
    if not groups_file.exists() or not membership_file.exists():
        LOGGER.warning("groups-file %s or membership-file %s does not exist", groups_file, membership_file)
        return
        
    try:
        membership_df = pd.read_csv(membership_file)
        if "membership" not in membership_df.columns or "sample" not in membership_df.columns:
            LOGGER.warning("membership file is missing columns 'membership' or 'sample'")
            return
            
        # Filter for cluster samples
        cluster_samples = membership_df.loc[membership_df["membership"] == cluster_id, "sample"]
        if cluster_samples.empty:
            LOGGER.warning("No samples found for cluster %d in membership file", cluster_id)
            return
            
        cluster_samples = cluster_samples.astype(str).tolist()
        
        # Load groups
        groups_df = pd.read_csv(groups_file)
        if "sample" not in groups_df.columns:
            LOGGER.warning("groups file is missing column 'sample'")
            return
            
        groups_df["sample"] = groups_df["sample"].astype(str)
        groups_subset = groups_df[groups_df["sample"].isin(cluster_samples)]
        if groups_subset.empty:
            LOGGER.warning("No matching groups data for cluster %d samples", cluster_id)
            return
            
        cell_type_col = "cell_type_adjusted" if "cell_type_adjusted" in groups_subset.columns else "cell_type"
        count_col = "count" if "count" in groups_subset.columns else ("0" if "0" in groups_subset.columns else None)
        if not count_col:
            LOGGER.warning("Could not find count column in groups file")
            return
            
        # Pivot to sample x cell_type matrix
        matrix_counts = (
            groups_subset.pivot_table(index="sample", columns=cell_type_col, values=count_col, aggfunc="sum")
            .fillna(0)
        )
        
        # Relative Abundance of (Counts + 1)
        counts_plus1 = matrix_counts + 1
        row_sums = counts_plus1.sum(axis=1).replace(0, np.nan)
        df_rel = counts_plus1.div(row_sums, axis=0).fillna(0)
        axis_label = "Rel. Abundance (Counts+1)"
        
        # Dataset map for coloring
        dataset_map = groups_subset[["sample", "dataset"]].drop_duplicates().set_index("sample")["dataset"]
        
        # Load network edges
        first_line = pd.read_csv(graph_file, sep="\t", nrows=1)
        if first_line.empty:
            LOGGER.info("No edges found in network file %s, skipping coabundance plots", graph_file)
            return
            
        has_header = "source" in first_line.columns
        if has_header:
            edges_df = pd.read_csv(graph_file, sep="\t")
        else:
            edges_df = pd.read_csv(graph_file, sep="\t", header=None)
            col_count = edges_df.shape[1]
            if col_count >= 2:
                cols = ["source", "target"]
                if col_count >= 3:
                    cols.append("MI")
                if col_count >= 4:
                    cols.append("pvalue")
                if col_count >= 5:
                    cols.append("pearson")
                if col_count >= 6:
                    cols.append("pearson_pvalue")
                if col_count >= 7:
                    cols.append("status")
                for i in range(len(cols), col_count):
                    cols.append(f"col_{i}")
                edges_df.columns = cols[:col_count]
            else:
                LOGGER.warning("Graph file has fewer than 2 columns, cannot identify edges")
                return
            
        if "source" in edges_df.columns:
            # Sort edges if MI or pearson exists
            if "MI" in edges_df.columns:
                edges_df = edges_df.sort_values(by="MI", ascending=False)
            elif "pearson" in edges_df.columns:
                edges_df = edges_df.sort_values(by="pearson", key=abs, ascending=False)
            
            # Select top 20 edges to plot
            top_edges = edges_df.head(20)
            
            for row in top_edges.itertuples(index=False):
                cell1 = str(row.source)
                cell2 = str(row.target)
                
                safe_cell1 = cell1.replace(" ", "_").replace("/", "-")
                safe_cell2 = cell2.replace(" ", "_").replace("/", "-")
                out_plot = output_dir / f"coabundance_{safe_cell1}_vs_{safe_cell2}.png"
                
                plot_pair_helper(df_rel, cell1, cell2, out_plot, dataset_map, axis_label)
        else:
            LOGGER.warning("Graph file does not contain 'source' column, cannot identify edges")
    except Exception as exc:
        LOGGER.error("Failed during coabundance plotting: %s", exc)

def plot_sample_pca(abundance_file: Path, membership_file: Path, output_dir: Path, time_label: str):
    LOGGER.info("Plotting sample PCA for %s using membership %s", abundance_file, membership_file)
    if not abundance_file.exists() or not membership_file.exists():
        LOGGER.warning("abundance-file %s or membership-file %s does not exist", abundance_file, membership_file)
        return
        
    try:
        # Load abundance matrix (sample x cell_type)
        df_ab = pd.read_csv(abundance_file, index_col=0)
        
        # Load membership
        df_mem = pd.read_csv(membership_file).set_index("sample")
        
        # Clean abundance df to have only cell types as features
        meta_cols = ["dataset", "membership"]
        features_cols = [c for c in df_ab.columns if c not in meta_cols]
        
        X = df_ab[features_cols].values
        
        if X.shape[0] < 2:
            LOGGER.warning("Not enough samples (%d) to run PCA", X.shape[0])
            return
            
        # Ensure we align memberships with abundance rows
        df_ab["membership"] = df_ab.index.map(df_mem["membership"]).fillna(-1).astype(int)
        memberships = df_ab["membership"]
        
        # Run PCA
        from sklearn.decomposition import PCA
        n_components = min(2, X.shape[1])
        if n_components < 2:
            LOGGER.warning("Not enough features (%d) to run PCA", X.shape[1])
            return
            
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X)
        
        # Plot
        fig, ax = plt.subplots(figsize=(8, 6))
        
        unique_memberships = sorted(memberships.unique())
        cmap = plt.get_cmap("tab20", len(unique_memberships))
        color_map = {m: cmap(i) for i, m in enumerate(unique_memberships)}
        colors = memberships.map(color_map)
        
        # Scatter plot
        ax.scatter(X_pca[:, 0], X_pca[:, 1], c=colors, alpha=0.8, s=40, edgecolor='k', linewidth=0.3)
        
        # Add labels and title
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title(f"Sample Abundance PCA (Painted by Cluster) - {time_label}")
        ax.grid(True, which="both", ls="--", alpha=0.2)
        
        # Add legend
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[m], label=f"Cluster {m}" if m != -1 else "Unknown", markersize=10) for m in unique_memberships]
        ax.legend(handles=handles, title="Ecotype Cluster", loc="best")
        
        pca_file = output_dir / f"sample_pca_abundance_{time_label}.png"
        fig.savefig(pca_file, dpi=200, bbox_inches='tight')
        plt.close(fig)
        LOGGER.info("Saved PCA abundance plot to %s", pca_file)
        
    except Exception as exc:
        LOGGER.error("Failed to generate PCA abundance plot: %s", exc)

def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Degree distribution plots (linear and semi-log)
    plot_degree_distribution(args.graph_file, args.output_dir, args.time_label)
    
    # 2. Coabundance plots (if groups and membership files are provided)
    if args.groups_file and args.membership_file:
        plot_cell_type_coabundance(args.graph_file, args.groups_file, args.membership_file, args.output_dir, args.time_label)
        
    # 3. PCA plots (if abundance and global-membership files are provided)
    if args.abundance_file and args.global_membership_file:
        plot_sample_pca(args.abundance_file, args.global_membership_file, args.output_dir, args.time_label)

if __name__ == "__main__":
    main()
