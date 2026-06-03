import os
import sys
import argparse
import json
import warnings
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import yaml

# Monkey-patch older pandas types used in pickled PyTorch scvi-tools checkpoints
class NumericModule:
    Int64Index = pd.Index
    Float64Index = pd.Index
    UInt64Index = pd.Index
    
import sys
sys.modules['pandas.core.indexes.numeric'] = NumericModule()

warnings.filterwarnings("ignore")

def map_gene_symbols_to_ensembl(adata, output_csv_path):
    import urllib.request
    import urllib.parse
    import json
    import time
    import pandas as pd
    
    print("Mapping gene symbols to Ensembl via MyGene.info API...")
    symbols = list(adata.var.index.dropna().unique())
    chunk_size = 200
    all_results = []
    
    for i in range(0, len(symbols), chunk_size):
        chunk = symbols[i:i+chunk_size]
        data = urllib.parse.urlencode({'q': ','.join(chunk), 'scopes': 'symbol', 'fields': 'ensembl.gene', 'species': 'human'}).encode('utf-8')
        req = urllib.request.Request('https://mygene.info/v3/query', data=data)
        
        success = False
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=15) as response:
                    res = json.loads(response.read().decode('utf-8'))
                    all_results.extend(res)
                    success = True
                    break
            except Exception as e:
                print(f"Attempt {attempt + 1}/5 failed for chunk starting at {i}: {e}")
                time.sleep(2 ** attempt)
                
        if not success:
            print(f"Warning: Failed to fetch mapping for chunk starting at {i}. Retrying with sub-chunks of 10...")
            sub_chunk_size = 10
            for j in range(0, len(chunk), sub_chunk_size):
                sub_chunk = chunk[j:j+sub_chunk_size]
                sub_data = urllib.parse.urlencode({'q': ','.join(sub_chunk), 'scopes': 'symbol', 'fields': 'ensembl.gene', 'species': 'human'}).encode('utf-8')
                sub_req = urllib.request.Request('https://mygene.info/v3/query', data=sub_data)
                
                sub_success = False
                for attempt in range(3):
                    try:
                        with urllib.request.urlopen(sub_req, timeout=10) as response:
                            res = json.loads(response.read().decode('utf-8'))
                            all_results.extend(res)
                            sub_success = True
                            break
                    except Exception as e:
                        time.sleep(1)
                if not sub_success:
                    print(f"Warning: Failed sub-chunk starting at {i + j}")
            
    mapping = {}
    for item in all_results:
        symbol = item.get('query')
        if 'ensembl' in item:
            ensembl = item['ensembl']
            if isinstance(ensembl, list):
                ensembl_id = ensembl[0].get('gene')
            else:
                ensembl_id = ensembl.get('gene')
            if ensembl_id:
                mapping[symbol] = ensembl_id
                
    mapping_df = pd.DataFrame(list(mapping.items()), columns=['symbol', 'ensembl_gene_id'])
    mapping_df.to_csv(output_csv_path, index=False)
    
    mapping_df.set_index('symbol', inplace=True)
    adata.var['symbol'] = adata.var.index
    mapori = adata.var.join(mapping_df, on='symbol')
    
    # avoid dropping all if duplicated, just keep first
    mapori = mapori[~mapori.index.duplicated(keep='first')]
    adata.var = mapori
    return adata

def load_config(config_path):
    if config_path.endswith('.yaml') or config_path.endswith('.yml'):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    elif config_path.endswith('.json'):
        with open(config_path, 'r') as f:
            return json.load(f)
    else:
        raise ValueError("Config must be .yaml or .json")

def main():
    parser = argparse.ArgumentParser(description="Run SCANVI Surgery")
    parser.add_argument("--config", type=str, required=True, help="Path to run config")
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_dir = cfg["run_dir"]
    os.makedirs(run_dir, exist_ok=True)

    scvi.settings.seed = cfg.get("seed", 0)
    torch.set_float32_matmul_precision(cfg.get("float32_matmul_precision", "high"))
    
    dataset_name = cfg["dataset_name"]
    name_short = cfg.get("dataset_name_short", dataset_name)

    print("Loading target dataset...")
    target_adata = sc.read_h5ad(cfg["input_h5ad"])

    # Basic QC
    target_adata.var['mito'] = target_adata.var.index.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(target_adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)

    max_cells = cfg.get("max_cells", None)
    if max_cells and target_adata.shape[0] > max_cells:
        print(f"Subsampling to {max_cells} cells...")
        sc.pp.subsample(target_adata, n_obs=max_cells)

    target_adata.obs.index = f'{name_short}_' + target_adata.obs.index.astype(str)
    target_adata.obs_names_make_unique()

    if cfg.get("map_ensembl", True):
        output_csv_path = os.path.join(run_dir, f'{name_short}_ensembl.csv')
        if not os.path.exists(output_csv_path):
            target_adata = map_gene_symbols_to_ensembl(target_adata, output_csv_path)
        else:
            print("Loading existing gene mapping CSV...")
            ensembl = pd.read_csv(output_csv_path)
            mapping_dict = dict(zip(ensembl['symbol'].astype(str), ensembl['ensembl_gene_id'].astype(str)))
            target_adata.var['symbol'] = target_adata.var.index.astype(str)
            target_adata.var['ensembl_gene_id'] = target_adata.var['symbol'].map(mapping_dict)
                
        # Drop genes without ensembl ID
        if 'ensembl_gene_id' in target_adata.var:
            target_adata = target_adata[:, target_adata.var['ensembl_gene_id'].notnull()]
            target_adata.var.index = target_adata.var['ensembl_gene_id'].astype(str)
            target_adata.var_names_make_unique()

    ref_path = cfg["reference_model"]
    
    if os.path.exists(os.path.join(ref_path, "luca_99perc_reference.h5ad")):
        print("Checking reference adata for overlap...")
        try:
            ref_adata = sc.read_h5ad(os.path.join(ref_path, "luca_99perc_reference.h5ad"), backed='r')
            intersect = np.intersect1d(target_adata.var.index, ref_adata.var.index)
            print(f"Number of overlapping genes: {len(intersect)}")
        except Exception as e:
            print(f"Could not read reference for overlap check: {e}")

    print("Preparing query anndata...")
    adata_query = scvi.model.SCANVI.prepare_query_anndata(
        adata=target_adata,
        reference_model=ref_path,
        inplace=False
    )
    
    # Ensure sample and dataset annotations
    if 'dataset' not in adata_query.obs:
        adata_query.obs['dataset'] = name_short
    if 'sample' not in adata_query.obs:
        if 'orig.ident' in adata_query.obs:
            adata_query.obs['sample'] = adata_query.obs['orig.ident']
        else:
            adata_query.obs['sample'] = adata_query.obs['dataset']
        
    print("Loading query data into SCANVI...")
    model = scvi.model.SCANVI.load_query_data(
        adata_query,
        ref_path,
        freeze_dropout=True,
    )
    
    # Setup training kwargs
    plan_kwargs = cfg.get("scanvi_plan_kwargs", {})
    train_kwargs = cfg.get("scanvi_train_params", {"max_epochs": 100})
    
    early_stopping = cfg.get("scanvi_early_stopping", False)
    if early_stopping:
        train_kwargs["early_stopping"] = True
        train_kwargs["early_stopping_monitor"] = "elbo_validation"
        train_kwargs["early_stopping_patience"] = 10
        train_kwargs["early_stopping_min_delta"] = 0.0

    print("Training SCANVI model on query...")
    model.train(
        plan_kwargs=plan_kwargs,
        **train_kwargs
    )
    
    print("Saving model...")
    model.save(os.path.join(run_dir, "scanvi_model"), overwrite=True)
    
    print("Getting latent representation...")
    adata_query.obsm["X_scVI"] = model.get_latent_representation()
    
    print("Predicting cell types...")
    adata_query.obs["predicted_cell_type"] = model.predict()
    
    print("Computing neighbors and UMAP...")
    sc.pp.neighbors(adata_query, use_rep="X_scVI")
    sc.tl.umap(adata_query)
    sc.tl.leiden(adata_query, flavor="igraph")
    
    print("Saving artifacts...")
    from pandas.api.types import is_string_dtype

    # Convert indexes and Arrow string arrays to standard numpy/pandas formats to avoid IORegistryError
    if hasattr(adata_query.obs, 'index') and adata_query.obs.index is not None:
        if 'Arrow' in str(type(adata_query.obs.index.array if hasattr(adata_query.obs.index, 'array') else adata_query.obs.index.values)) or 'Arrow' in str(adata_query.obs.index.dtype):
            adata_query.obs.index = pd.Index(adata_query.obs.index.tolist(), dtype=object)
            
    if hasattr(adata_query.var, 'index') and adata_query.var.index is not None:
        if 'Arrow' in str(type(adata_query.var.index.array if hasattr(adata_query.var.index, 'array') else adata_query.var.index.values)) or 'Arrow' in str(adata_query.var.index.dtype):
            adata_query.var.index = pd.Index(adata_query.var.index.tolist(), dtype=object)
        
    for col in adata_query.obs.columns:
        # Pre-emptively convert string/object columns to category
        if is_string_dtype(adata_query.obs[col]):
            adata_query.obs[col] = adata_query.obs[col].astype('category')
            
        if isinstance(adata_query.obs[col].dtype, pd.CategoricalDtype):
            cats = adata_query.obs[col].cat.categories
            if 'Arrow' in str(type(cats.array if hasattr(cats, 'array') else cats.values)) or 'Arrow' in str(cats.dtype):
                adata_query.obs[col] = adata_query.obs[col].cat.rename_categories(
                    pd.Index(cats.tolist(), dtype=object)
                )
            
    for col in adata_query.var.columns:
        # Pre-emptively convert string/object columns to category
        if is_string_dtype(adata_query.var[col]):
            adata_query.var[col] = adata_query.var[col].astype('category')
            
        if isinstance(adata_query.var[col].dtype, pd.CategoricalDtype):
            cats = adata_query.var[col].cat.categories
            if 'Arrow' in str(type(cats.array if hasattr(cats, 'array') else cats.values)) or 'Arrow' in str(cats.dtype):
                adata_query.var[col] = adata_query.var[col].cat.rename_categories(
                    pd.Index(cats.tolist(), dtype=object)
                )

    adata_query.write_h5ad(os.path.join(run_dir, "query_latent.h5ad"))
    
    # Output CSV mapping original cell types to predicted cell types if original exists
    cell_type_key = cfg.get("cell_type_key", "cell_type_tumor")
    if cell_type_key in adata_query.obs:
        df = adata_query.obs[[cell_type_key, "predicted_cell_type", "leiden", "sample"]].copy()
    else:
        df = adata_query.obs[["predicted_cell_type", "leiden", "sample"]].copy()
        
    df.to_csv(os.path.join(run_dir, f"{name_short}_predicted.csv"))
    
    print("Surgery complete.")
    
if __name__ == "__main__":
    raise SystemExit(main())
