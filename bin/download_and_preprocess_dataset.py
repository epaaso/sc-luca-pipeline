#!/usr/bin/env python3
import os
import sys
import argparse
import urllib.request
import pathlib
import subprocess
import shutil
import re
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

def download_file(url, local_path):
    local_path = pathlib.Path(local_path)
    if local_path.exists():
        print(f"File already exists: {local_path}")
        return local_path
    
    local_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} to {local_path}...")
    
    # Simple block-based download with progress print
    temp_path = local_path.with_suffix(".tmp")
    try:
        urllib.request.urlretrieve(url, temp_path)
        temp_path.rename(local_path)
        print(f"Downloaded {local_path} successfully.")
    except Exception as e:
        if temp_path.exists():
            temp_path.unlink()
        print(f"Error downloading {url}: {e}")
        raise e
    return local_path

def download_parallel(urls_paths, cpus=4):
    with ThreadPoolExecutor(max_workers=cpus) as executor:
        futures = [executor.submit(download_file, url, path) for url, path in urls_paths]
        for future in futures:
            future.result()

def run_r_script(code_str):
    print("Running R script for data preparation/conversion...")
    result = subprocess.run(["Rscript", "-e", code_str], capture_output=True, text=True)
    if result.returncode != 0:
        print("R script stdout:")
        print(result.stdout)
        print("R script stderr:")
        print(result.stderr)
        raise RuntimeError(f"R script failed with exit code {result.returncode}")
    print("R script completed successfully.")

def map_tnm_to_roman(tnm_stage):
    if not isinstance(tnm_stage, str):
        return "Unknown stage"
    stage_mappings = {
        r"^TisN0M0$": "0",
        r"^T1[abc]?N0M0$": "I",
        r"^T2[abc]?N0M0$": "II",
        r"^T2[abc]?N0/1$": "II",
        r"^T1[abc]?N1M0$": "II",
        r"^T2[abc]?N1M0$": "II",
        r"^T3N[01]M0$": "II",
        r"^T3M0N[01]$": "II",
        r"^T1[abc]?N2M0$": "III",
        r"^T2[abc]?N2M0$": "III",
        r"^T3N2M0$": "III",
        r"^T4N[012]M0$": "III",
        r"^T4N0$": "III",
        r"^T[1-4][abc]?N3M0$": "III",
        r"^T[1-4][abc]?N[0-3]M1[abc]$": "IV",
    }
    for pattern, stage in stage_mappings.items():
        if re.match(pattern, tnm_stage):
            return stage
    return "Unknown stage"

def concatenate_datasets(directory, sample_ids=None):
    adatas = {}
    required_files = {"barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"}
    
    directory = pathlib.Path(directory)
    if sample_ids is None:
        subdirectories = sorted(d.name for d in directory.iterdir() if d.is_dir())
    else:
        subdirectories = list(sample_ids)
        
    for subdirectory in subdirectories:
        subdirectory_path = directory / subdirectory
        if not subdirectory_path.is_dir():
            raise FileNotFoundError(f"Missing sample directory: {subdirectory_path}")
            
        files = set(f.name for f in subdirectory_path.iterdir())
        if not required_files.issubset(files):
            nested_path = subdirectory_path / "filtered_feature_bc_matrix"
            nested_files = set(f.name for f in nested_path.iterdir()) if nested_path.is_dir() else set()
            if required_files.issubset(nested_files):
                subdirectory_path = nested_path
            else:
                missing = sorted(required_files - files)
                raise FileNotFoundError(f"Missing 10x files in {subdirectory_path}: {missing}")
                
        adata = sc.read_10x_mtx(subdirectory_path, var_names="gene_symbols", make_unique=True)
        adata.obs_names = [f"{subdirectory}_{barcode}" for barcode in adata.obs_names]
        adatas[subdirectory] = adata
        
    if not adatas:
        raise ValueError(f"No 10x sample directories found in {directory}")
        
    adata = ad.concat(adatas, label="sample", index_unique=None)
    return adata

def bishoff_pipeline(data_dir, output_file):
    print("Executing Bishoff dataset preparation...")
    bishoff_dir = data_dir / "Bishoff2021"
    bishoff_dir.mkdir(parents=True, exist_ok=True)
    
    bishoff_samples = ["p019t", "p030t", "p033t"]
    tenx_filenames = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]
    codeocean_tree_base = "https://codeocean.com/capsule/8321305/tree/v1/data"
    codeocean_download_base = "https://files.codeocean.com/files/verified/fcd176f7-ea81-4e4c-a556-86f0683abb47_v1.0/data"
    
    urls_paths = []
    for sample in bishoff_samples:
        sample_dir = bishoff_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        for filename in tenx_filenames:
            local_path = sample_dir / filename
            tree_url = f"{codeocean_tree_base}/cellranger/{sample}/filtered_feature_bc_matrix/{filename}"
            download_url = tree_url.replace(codeocean_tree_base, codeocean_download_base) + "?download="
            urls_paths.append((download_url, local_path))
            
    download_parallel(urls_paths)
    
    metadata_path = bishoff_dir / "patients_metadata.xlsx"
    metadata_url = f"{codeocean_download_base}/metadata/patients_metadata.xlsx?download="
    download_file(metadata_url, metadata_path)
    
    print("Concatenating Bishoff samples...")
    adata = concatenate_datasets(bishoff_dir, sample_ids=bishoff_samples)
    
    print("Merging Bishoff metadata...")
    meta = pd.read_excel(metadata_path).rename(columns={"sample_id": "sample"})
    adata.obs["barcode"] = adata.obs.index
    meta_cols = [c for c in meta.columns if c != "sample"]
    adata.obs = adata.obs.drop(columns=[c for c in meta_cols if c in adata.obs.columns], errors="ignore")
    adata.obs = adata.obs.merge(meta, how="left", on="sample").set_index("barcode")
    adata.obs.index.name = None
    
    datetime_cols = adata.obs.select_dtypes(include=["datetime64[ns]", "datetimetz"]).columns
    for col in datetime_cols:
        adata.obs[col] = adata.obs[col].dt.strftime("%Y-%m-%d %H:%M:%S")
        
    for col in adata.obs.select_dtypes(include=["category"]).columns:
        categories = adata.obs[col].cat.categories
        if pd.api.types.is_datetime64_any_dtype(categories):
            adata.obs[col] = adata.obs[col].astype("object")
            adata.obs[col] = adata.obs[col].where(adata.obs[col].isna(), adata.obs[col].astype(str))
            
    ad.settings.allow_write_nullable_strings = True
    
    print("Running Quality Control & Filtering...")
    # min_counts: 1000, max_counts: 100000, min_genes: 500, max_genes: 10000, max_pct_mito: 20
    adata.var["mito"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    sc.pp.filter_genes(adata, min_counts=2)
    sc.pp.filter_cells(adata, min_counts=1000)
    sc.pp.filter_cells(adata, max_counts=100000)
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_cells(adata, max_genes=10000)
    adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
    
    adata.write_h5ad(output_file)
    print(f"Bishoff dataset successfully prepared: {output_file}")

def hu_pipeline(data_dir, output_file):
    print("Executing Hu dataset preparation...")
    hu_dir = data_dir / "Hu2023"
    hu_dir.mkdir(parents=True, exist_ok=True)
    
    matrix_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207422&format=file&file=GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    metadata_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207422&format=file&file=GSE207422%5FNSCLC%5FscRNAseq%5Fmetadata%2Exlsx"
    
    matrix_gz = hu_dir / "umis.txt.gz"
    matrix_txt = hu_dir / "umis.txt"
    metadata_path = hu_dir / "metadata.xlsx"
    
    download_file(matrix_url, matrix_gz)
    download_file(metadata_url, metadata_path)
    
    if not matrix_txt.exists():
        print("Extracting matrix UMI counts...")
        import gzip
        with gzip.open(matrix_gz, "rb") as f_in:
            with open(matrix_txt, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
                
    # Direct Transpose using pandas (if memory permits, chunked transposition/tsv reader is used in get_data)
    print("Reading and transposing UMI matrix (this may take a few minutes)...")
    df = pd.read_csv(matrix_txt, sep="\t", index_col=0)
    df_t = df.T
    adata = ad.AnnData(df_t)
    
    print("Merging Hu metadata...")
    meta = pd.read_excel(metadata_path)
    adata.obs = adata.obs.merge(meta, left_index=True, right_on="Cell_Index", how="left").set_index("Cell_Index")
    adata.obs.index.name = None
    
    print("Running Quality Control & Filtering...")
    # min_genes: 500, max_pct_mito: 20
    adata.var["mito"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    sc.pp.filter_genes(adata, min_counts=2)
    sc.pp.filter_cells(adata, min_genes=500)
    adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
    
    adata.write_h5ad(output_file)
    print(f"Hu dataset successfully prepared: {output_file}")

def zuani_pipeline(data_dir, output_file):
    print("Executing Zuani dataset preparation...")
    zuani_dir = data_dir / "Zuani2024"
    zuani_dir.mkdir(parents=True, exist_ok=True)
    
    sdrf_url = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/E-MTAB-13526.sdrf.txt"
    sdrf_path = zuani_dir / "E-MTAB-13526.sdrf.txt"
    download_file(sdrf_url, sdrf_path)
    
    zuani_samples = pd.read_csv(sdrf_path, sep="\t")
    facs_filter = ["sorted tumor cells", "sorted epithelial cells", "unsorted lung cells"]
    zuani_samples = zuani_samples[zuani_samples["Factor Value[FACS]"].isin(facs_filter)]
    zuani_samples = zuani_samples[zuani_samples["Factor Value[sampling site]"] == "tumor"]
    
    zuani_patients = zuani_samples["Source Name"].unique()
    base_url = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files"
    
    urls_paths = []
    for pat in zuani_patients:
        # Download barcodes, features, matrix for each patient sample
        for ext in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
            filename = f"{pat}-{ext.replace('.tsv.gz', '.tsv').replace('.mtx.gz', '.mtx')}.gz"
            # Map EBI biostudies format to standard cellranger filenames
            # get_data.ipynb expects them to be renamed to barcodes.tsv.gz, etc inside directories
            patient_dir = zuani_dir / pat
            patient_dir.mkdir(parents=True, exist_ok=True)
            local_filename = ext
            urls_paths.append((f"{base_url}/{filename}", patient_dir / local_filename))
            
    print(f"Downloading raw data files for {len(zuani_patients)} Zuani patients...")
    download_parallel(urls_paths)
    
    print("Concatenating Zuani datasets...")
    adata = concatenate_datasets(zuani_dir, sample_ids=zuani_patients)
    
    print("Mapping Zuani metadata...")
    # Normalize sample names to merge metadata
    adata.obs["sample"] = adata.obs["sample"].str[:-1]
    zuani_samples_unique = zuani_samples.drop_duplicates(subset=["Source Name"], keep="first")
    adata.obs = pd.merge(adata.obs, zuani_samples_unique.loc[:, ["Source Name", "Characteristics[tumor grading]"]],
                         left_on="sample", right_on="Source Name", how="left")
    
    adata = adata[~adata.obs["Characteristics[tumor grading]"].isna()]
    adata.obs["stage"] = adata.obs["Characteristics[tumor grading]"].apply(map_tnm_to_roman)
    
    print("Running Quality Control & Filtering...")
    # min_counts: 400, max_counts: 100000, min_genes: 180, max_genes: 6000, max_pct_mito: 20
    adata.var["mito"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    sc.pp.filter_genes(adata, min_counts=2)
    sc.pp.filter_cells(adata, min_counts=400)
    sc.pp.filter_cells(adata, max_counts=100000)
    sc.pp.filter_cells(adata, min_genes=180)
    sc.pp.filter_cells(adata, max_genes=6000)
    adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
    
    adata.write_h5ad(output_file)
    print(f"Zuani dataset successfully prepared: {output_file}")

def deng_pipeline(data_dir, output_file):
    print("Executing Deng dataset preparation...")
    deng_dir = data_dir / "Deng2024"
    deng_dir.mkdir(parents=True, exist_ok=True)
    
    rds_url = "https://api.figshare.com/v2/file/download/46617667"
    rds_path = deng_dir / "dengEcoTyper.rds"
    download_file(rds_url, rds_path)
    
    h5seurat_path = deng_dir / "deng.h5Seurat"
    h5ad_raw_path = deng_dir / "dengRNA.h5ad"
    
    r_code = f"""
    library(Seurat)
    library(SeuratDisk)
    deng <- readRDS("{rds_path.resolve()}")
    deng = UpdateSeuratObject(deng)
    cells_to_exclude <- Cells(deng)[grepl("P05|N", deng@meta.data$orig.ident)]
    deng <- subset(deng, cells = setdiff(Cells(deng), cells_to_exclude))
    SaveH5Seurat(deng, filename="{h5seurat_path.resolve()}", overwrite=TRUE)
    Convert("{h5seurat_path.resolve()}", dest = "{h5ad_raw_path.resolve()}", assay = "RNA", overwrite=TRUE)
    """
    run_r_script(r_code)
    
    print("Reading converted Deng AnnData...")
    adata = sc.read_h5ad(h5ad_raw_path)
    
    print("Running Quality Control & Filtering...")
    # min_genes: 200, max_genes: 6000, max_pct_mito: 20
    adata.var["mito"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
    sc.pp.filter_genes(adata, min_counts=2)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, max_genes=6000)
    adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
    
    adata.write_h5ad(output_file)
    print(f"Deng dataset successfully prepared: {output_file}")

def extended_atlas_pipeline(data_dir, output_file):
    print("Executing Extended Atlas dataset download...")
    extended_url = "https://datasets.cellxgene.cziscience.com/173984ce-d33a-46b6-ae96-4be47f6c67e8.h5ad"
    download_file(extended_url, output_file)
    print(f"Extended Atlas successfully downloaded: {output_file}")

def extended_tumor_hvg_pipeline(data_dir, output_file):
    print("Executing Extended Tumor HVG Atlas dataset preparation...")
    extended_path = data_dir / "extended.h5ad"
    extended_url = "https://datasets.cellxgene.cziscience.com/173984ce-d33a-46b6-ae96-4be47f6c67e8.h5ad"
    download_file(extended_url, extended_path)
    
    print("Loading full extended atlas...")
    adata = sc.read_h5ad(extended_path)
    
    print("Subsetting to tumor_primary...")
    adata = adata[adata.obs.origin == 'tumor_primary'].copy()
    
    print("Filtering uicc_stage...")
    stages = ['I', 'II', 'III', 'III or IV', 'IV']
    adata = adata[adata.obs.uicc_stage.isin(stages)].copy()
    
    print("Filtering study...")
    studies = ['Goveia_Carmeliet_2020', 'Leader_Merad_2021', 'Guo_Zhang_2018']
    adata = adata[~adata.obs.study.isin(studies)].copy()
    
    print("Filtering highly variable genes...")
    hvg_mask = (adata.var.is_highly_variable == True) | (adata.var.is_highly_variable == 'True')
    adata = adata[:, hvg_mask].copy()
    
    if hasattr(adata, 'raw') and adata.raw is not None:
        del adata.raw
        
    if 'counts_length_scaled' in adata.layers:
        del adata.layers['counts_length_scaled']
        
    if 'count' in adata.layers:
        adata.X = adata.layers['count'].copy()
        
    print("Writing subsetted atlas h5ad...")
    adata.write_h5ad(output_file)
    print(f"Extended Tumor HVG Atlas successfully prepared: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Download and preprocess datasets automatically.")
    parser.add_argument("--dataset", required=True, choices=["Bishoff", "Hu", "Zuani", "Deng", "extended", "extended_tumor_hvg"], help="Dataset to prepare.")
    parser.add_argument("--output", required=True, help="Output path for the prepared .h5ad file.")
    parser.add_argument("--data-dir", default="./data_downloads", help="Directory for temporary raw downloads.")
    args = parser.parse_args()
    
    output_path = pathlib.Path(args.output)
    if output_path.exists():
        print(f"Output file already exists at {output_path}. Skipping preparation.")
        return 0
        
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data_dir = pathlib.Path(args.data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    
    if args.dataset == "Bishoff":
        bishoff_pipeline(data_dir, output_path)
    elif args.dataset == "Hu":
        hu_pipeline(data_dir, output_path)
    elif args.dataset == "Zuani":
        zuani_pipeline(data_dir, output_path)
    elif args.dataset == "Deng":
        deng_pipeline(data_dir, output_path)
    elif args.dataset == "extended":
        extended_atlas_pipeline(data_dir, output_path)
    elif args.dataset == "extended_tumor_hvg":
        extended_tumor_hvg_pipeline(data_dir, output_path)
        
    return 0

if __name__ == "__main__":
    sys.exit(main())
