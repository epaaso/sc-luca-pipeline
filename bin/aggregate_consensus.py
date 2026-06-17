#!/usr/bin/env python3
"""Aggregate marker genes and Hallmark GSEA scores across datasets for each ecotype."""

import argparse
import glob
import logging
import os
import re
import sys
from typing import List, Dict, Tuple
import numpy as np
try:
    import numpy.rec
except ImportError:
    try:
        import numpy.core.records as records
        sys.modules['numpy.rec'] = records
    except ImportError:
        pass
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate markers and GSEA scores for every ecotype.")
    parser.add_argument("--marker-dir", required=True, help="Directory containing npy marker files")
    parser.add_argument("--gsea-dir", required=True, help="Directory containing GSEA csv files")
    parser.add_argument("--output-dir", required=True, help="Output directory for consensus results")
    return parser.parse_args()

def aggregate_markers(stage: str, cluster: int, marker_files: List[str], output_dir: str) -> None:
    logging.info("Aggregating marker genes for ecotype %s cluster %02d (from %d files)", stage, cluster, len(marker_files))
    
    # cell_type -> {gene -> list of scores}
    all_cell_types: Dict[str, Dict[str, List[float]]] = {}
    
    for f in marker_files:
        try:
            d = np.load(f, allow_pickle=True).item()
            if "names" not in d or "scores" not in d:
                logging.warning("File %s does not contain 'names' or 'scores' keys, skipping", f)
                continue
            
            # Get cell types
            cell_types = list(d["scores"].dtype.names)
            for ct in cell_types:
                if ct not in all_cell_types:
                    all_cell_types[ct] = {}
                
                names_arr = d["names"][ct]
                scores_arr = d["scores"][ct]
                for gene, score in zip(names_arr, scores_arr):
                    if gene not in all_cell_types[ct]:
                        all_cell_types[ct][gene] = []
                    all_cell_types[ct][gene].append(float(score))
        except Exception as exc:
            logging.error("Error loading marker file %s: %s", f, exc)
            
    if not all_cell_types:
        logging.warning("No cell types found to aggregate for %s cluster %02d", stage, cluster)
        return

    # Calculate average scores and sort
    consensus_data: Dict[str, List[Tuple[str, float]]] = {}
    for ct, gene_scores in all_cell_types.items():
        avg_scores = []
        for gene, scores in gene_scores.items():
            avg_scores.append((gene, float(np.mean(scores))))
        # Sort by average AUC in descending order
        avg_scores.sort(key=lambda x: x[1], reverse=True)
        consensus_data[ct] = avg_scores

    # 1. Save top 10 genes and scores as CSV
    top10_rows = []
    for ct, genes in consensus_data.items():
        for rank, (gene, score) in enumerate(genes[:10], 1):
            top10_rows.append({
                "cell_type": ct,
                "rank": rank,
                "gene": gene,
                "mean_auc": score
            })
    
    top10_df = pd.DataFrame(top10_rows)
    csv_path = os.path.join(output_dir, f"consensus_markers_{stage}_cluster_{cluster}_top10.csv")
    top10_df.to_csv(csv_path, index=False)
    logging.info("Saved top 10 consensus markers to %s", csv_path)

    # 2. Save structured numpy file for compatibility
    max_len = max(len(genes) for genes in consensus_data.values())
    
    # We pad shorter lists if needed (though they should be same length)
    cell_types = list(consensus_data.keys())
    dtypes_names = [(ct, "O") for ct in cell_types]
    dtypes_scores = [(ct, "<f8") for ct in cell_types]
    
    names_list = []
    scores_list = []
    for ct in cell_types:
        genes = consensus_data[ct]
        ct_names = [g[0] for g in genes] + [""] * (max_len - len(genes))
        ct_scores = [g[1] for g in genes] + [0.5] * (max_len - len(genes))
        names_list.append(ct_names)
        scores_list.append(ct_scores)
        
    names_arr = np.array(list(zip(*names_list)), dtype=dtypes_names)
    scores_arr = np.array(list(zip(*scores_list)), dtype=dtypes_scores)
    
    npy_dict = {
        "params": {"reference": "tumorall"},
        "names": names_arr,
        "scores": scores_arr
    }
    npy_path = os.path.join(output_dir, f"consensus_markers_{stage}_cluster_{cluster}.npy")
    np.save(npy_path, npy_dict, allow_pickle=True)
    logging.info("Saved consensus markers numpy object to %s", npy_path)

    # 3. Plot top 10 consensus markers
    num_types = len(cell_types)
    rows = (num_types + 1) // 2
    fig, axs = plt.subplots(rows, 2, figsize=(16, 4.5 * rows))
    if num_types == 1:
        axs = np.array([axs, axs]) # make it 1D
    axs = axs.ravel()
    
    for i, ct in enumerate(cell_types):
        top_genes = consensus_data[ct][:10]
        genes_names = [g[0] for g in top_genes]
        genes_scores = [g[1] for g in top_genes]
        
        # Plot horizontal bar chart (highest score at top)
        y_pos = np.arange(len(genes_names))
        axs[i].barh(y_pos, genes_scores, color=sns.color_palette("viridis", 10), edgecolor="k", linewidth=0.5)
        axs[i].set_yticks(y_pos)
        axs[i].set_yticklabels(genes_names, fontsize=8)
        axs[i].invert_yaxis()  # top-down ranking
        axs[i].set_xlim(0.4, 1.0)
        axs[i].set_xlabel("Mean AUC", fontsize=8)
        axs[i].set_title(f"{ct}", fontsize=10, fontweight="bold")
        axs[i].grid(True, axis="x", ls="--", alpha=0.3)
        
    # Delete unused subplots
    if len(axs) > num_types:
        for j in range(num_types, len(axs)):
            fig.delaxes(axs[j])
            
    plt.suptitle(f"Consensus Top 10 Marker Genes - Ecotype {stage} Cluster {cluster}", fontsize=14, fontweight="bold", y=0.99)
    plt.tight_layout()
    plot_path = os.path.join(output_dir, f"consensus_markers_{stage}_cluster_{cluster}.png")
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logging.info("Saved consensus markers plot to %s", plot_path)

def aggregate_gsea(stage: str, cluster: int, gsea_dir: str, output_dir: str) -> None:
    pattern = os.path.join(gsea_dir, "**", f"heatmap_{stage}_*_cluster_{cluster}_auc_tumorall_{stage}.csv")
    gsea_files = glob.glob(pattern, recursive=True)
    
    if not gsea_files:
        logging.warning("No GSEA files found for ecotype %s cluster %02d, skipping GSEA aggregation", stage, cluster)
        return
        
    logging.info("Aggregating Hallmark GSEA scores for ecotype %s cluster %02d (from %d files)", stage, cluster, len(gsea_files))
    
    dfs = []
    for f in gsea_files:
        try:
            df = pd.read_csv(f, index_col=0)
            if df.empty:
                continue
            dfs.append(df)
        except Exception as exc:
            logging.error("Error loading GSEA file %s: %s", f, exc)
            
    if not dfs:
        logging.warning("No valid GSEA dataframes loaded for %s cluster %02d", stage, cluster)
        return
        
    # Concatenate and compute mean of NES scores grouped by cell type index
    combined_df = pd.concat(dfs)
    consensus_df = combined_df.groupby(combined_df.index).mean()
    
    # 1. Save GSEA table
    csv_path = os.path.join(output_dir, f"consensus_gsea_{stage}_cluster_{cluster}.csv")
    consensus_df.to_csv(csv_path)
    logging.info("Saved consensus GSEA table to %s", csv_path)
    
    # 2. Plot GSEA heatmap
    plt.figure(figsize=(16, max(6, 0.4 * len(consensus_df))))
    sns.heatmap(
        consensus_df, 
        cmap="coolwarm", 
        center=0,
        cbar_kws={"label": "Mean NES"},
        linewidths=0.5,
        linecolor="gray",
        xticklabels=True,
        yticklabels=True
    )
    plt.title(f"Consensus Hallmark GSEA Mean NES - Ecotype {stage} Cluster {cluster}", fontsize=14, fontweight="bold")
    plt.xlabel("Hallmark Pathways", fontsize=10)
    plt.ylabel("Cell Types", fontsize=10)
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=8)
    
    plot_path = os.path.join(output_dir, f"consensus_gsea_{stage}_cluster_{cluster}.png")
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close()
    logging.info("Saved consensus GSEA heatmap plot to %s", plot_path)

def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    stages = ["I-II", "III-IV"]
    clusters = [0, 1, 2]
    
    for stage in stages:
        for cluster in clusters:
            # Find all marker npy files for this stage and cluster recursively under marker_dir
            pattern = os.path.join(args.marker_dir, "**", f"{stage}_*_cluster_{cluster}_auc.npy")
            marker_files = glob.glob(pattern, recursive=True)
                
            if marker_files:
                aggregate_markers(stage, cluster, marker_files, args.output_dir)
                aggregate_gsea(stage, cluster, args.gsea_dir, args.output_dir)
            else:
                logging.info("No marker files found for stage %s cluster %02d, skipping", stage, cluster)
                
    logging.info("Consensus aggregation completed successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
