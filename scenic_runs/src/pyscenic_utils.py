"""
pySCENIC Utilities for Organoid Analysis
Adapted from the original pyscenic_celline/ps_utils.py

This module contains utility functions for:
- Subsampling strategies
- HVG selection
- Network inference with Dask
"""

import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import random
import socket
from arboreto.algo import grnboost2
from distributed import LocalCluster, Client
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df
from pyscenic.aucell import aucell
import os
import warnings
warnings.filterwarnings('ignore')


def subsample_adata_strat(adata, subsample_column, subsample_ncells, strat=None, seed=42):
    """
    Subsample cells from AnnData object with stratification option.
    
    Args:
        adata: AnnData object
        subsample_column: Column to use for stratification
        subsample_ncells: Number of cells to subsample
        strat: Stratification strategy (None for random)
        seed: Random seed
        
    Returns:
        Subsampled AnnData object
    """
    np.random.seed(seed)
    random.seed(seed)
    
    if strat is None:
        # Random subsampling
        if subsample_ncells >= adata.n_obs:
            return adata.copy()
        
        indices = np.random.choice(adata.n_obs, subsample_ncells, replace=False)
        return adata[indices].copy()
    
    else:
        # Stratified subsampling
        if subsample_column not in adata.obs.columns:
            raise ValueError(f"Column {subsample_column} not found in adata.obs")
        
        # Get unique values and their counts
        value_counts = adata.obs[subsample_column].value_counts()
        
        # Calculate proportional sampling
        total_cells = len(adata.obs)
        selected_indices = []
        
        for value in value_counts.index:
            value_mask = adata.obs[subsample_column] == value
            value_indices = np.where(value_mask)[0]
            
            # Calculate how many cells to sample from this group
            proportion = len(value_indices) / total_cells
            n_sample = int(proportion * subsample_ncells)
            
            # Ensure we don't sample more than available
            n_sample = min(n_sample, len(value_indices))
            
            if n_sample > 0:
                sampled = np.random.choice(value_indices, n_sample, replace=False)
                selected_indices.extend(sampled)
        
        # If we haven't reached the target, randomly sample remaining
        if len(selected_indices) < subsample_ncells:
            remaining_indices = set(range(adata.n_obs)) - set(selected_indices)
            remaining_needed = subsample_ncells - len(selected_indices)
            additional = np.random.choice(list(remaining_indices), 
                                        min(remaining_needed, len(remaining_indices)), 
                                        replace=False)
            selected_indices.extend(additional)
        
        return adata[selected_indices].copy()


def hvg_adata(adata, n_top_genes=3000, flavor='seurat_v3'):
    """
    Select highly variable genes from AnnData object.
    
    Args:
        adata: AnnData object
        n_top_genes: Number of top HVGs to select
        flavor: Method for HVG selection
        
    Returns:
        AnnData object with HVGs
    """
    adata_hvg = adata.copy()
    
    # Calculate HVGs
    sc.pp.highly_variable_genes(adata_hvg, 
                               n_top_genes=n_top_genes, 
                               flavor=flavor)
    
    # Filter to HVGs
    adata_hvg = adata_hvg[:, adata_hvg.var.highly_variable].copy()
    
    return adata_hvg


def find_free_port():
    """
    Find a free port for Dask scheduler.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        s.listen(1)
        port = s.getsockname()[1]
    return port


def run_infer_partial_network(expression_data, tf_names, seed=42, n_workers=8, 
                             memory_limit="8e9", client=None):
    """
    Run GRNBoost2 network inference.
    
    Args:
        expression_data: Gene expression DataFrame
        tf_names: List of transcription factor names
        seed: Random seed
        n_workers: Number of workers
        memory_limit: Memory limit per worker
        client: Dask client (optional)
        
    Returns:
        Adjacency DataFrame
    """
    if client is None:
        # Create local cluster
        port = find_free_port()
        cluster = LocalCluster(
            n_workers=n_workers,
            threads_per_worker=2,
            memory_limit=memory_limit,
            scheduler_port=port,
            silence_logs=False
        )
        client = Client(cluster)
    
    try:
        # Run GRNBoost2
        adjacencies = grnboost2(
            expression_data=expression_data,
            tf_names=tf_names,
            client_or_address=client,
            seed=seed,
            verbose=True
        )
        
        return adjacencies
        
    finally:
        if 'cluster' in locals():
            client.close()
            cluster.close()


def run_pyscenic_workflow(adata, tf_names, ranking_dbs, motif_annotations,
                         output_dir, seed=42, hvg_genes=3000,
                         subsample_column=None, subsample_ncells=None, 
                         strat=None, n_workers=8, memory_limit="8e9"):
    """
    Complete pySCENIC workflow implementation.
    
    Args:
        adata: AnnData object with expression data
        tf_names: List of transcription factor names
        ranking_dbs: List of ranking database paths
        motif_annotations: Path to motif annotations
        output_dir: Output directory
        seed: Random seed
        hvg_genes: Number of HVGs to use
        subsample_column: Column for stratified subsampling
        subsample_ncells: Number of cells to subsample
        strat: Stratification strategy
        n_workers: Number of workers
        memory_limit: Memory limit per worker
        
    Returns:
        Tuple of (adjacencies, modules, regulons, auc_matrix)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Starting pySCENIC workflow with seed {seed}")
    print(f"Initial data shape: {adata.shape}")
    
    # Step 1: Subsampling (if requested)
    if subsample_ncells is not None and subsample_ncells < adata.n_obs:
        print(f"Subsampling to {subsample_ncells} cells...")
        adata = subsample_adata_strat(adata, subsample_column, subsample_ncells, strat, seed)
        print(f"After subsampling: {adata.shape}")
    
    # Step 2: HVG selection
    if hvg_genes is not None and hvg_genes < adata.n_vars:
        print(f"Selecting {hvg_genes} HVGs...")
        adata = hvg_adata(adata, n_top_genes=hvg_genes)
        print(f"After HVG selection: {adata.shape}")
    
    # Convert to DataFrame
    expression_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
        index=adata.obs.index,
        columns=adata.var.index
    )
    
    # Step 3: Network inference (GRNBoost2)
    print("Running GRNBoost2...")
    adjacencies = run_infer_partial_network(
        expression_df, tf_names, seed=seed, 
        n_workers=n_workers, memory_limit=memory_limit
    )
    
    # Save adjacencies
    adj_file = os.path.join(output_dir, f"adjacencies_seed{seed}.csv")
    adjacencies.to_csv(adj_file, sep='\t', index=False)
    print(f"Saved adjacencies: {adj_file}")
    
    # Step 4: Module discovery
    print("Discovering modules...")
    modules = list(modules_from_adjacencies(adjacencies, expression_df))
    
    # Save modules
    modules_file = os.path.join(output_dir, f"modules_seed{seed}.pkl")
    import pickle
    with open(modules_file, 'wb') as f:
        pickle.dump(modules, f)
    print(f"Saved modules: {modules_file}")
    
    # Step 5: Motif enrichment
    print("Running motif enrichment...")
    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    
    # Load ranking databases
    dbs = [RankingDatabase(fname=db) for db in ranking_dbs]
    
    # Run pruning
    regulons_df = prune2df(dbs, modules, motif_annotations, num_workers=n_workers)
    
    # Save regulons
    regulons_file = os.path.join(output_dir, f"regulons_seed{seed}.csv")
    regulons_df.to_csv(regulons_file, sep='\t', index=False)
    print(f"Saved regulons: {regulons_file}")
    
    # Convert to regulon objects
    from pyscenic.prune import df2regulons
    regulons = df2regulons(regulons_df)
    
    # Save regulon objects
    regulons_pkl = os.path.join(output_dir, f"regulons_seed{seed}.pkl")
    with open(regulons_pkl, 'wb') as f:
        pickle.dump(regulons, f)
    print(f"Saved regulon objects: {regulons_pkl}")
    
    # Step 6: AUCell scoring
    print("Calculating AUCell scores...")
    auc_matrix = aucell(expression_df, regulons, num_workers=n_workers)
    
    # Save AUCell matrix
    auc_file = os.path.join(output_dir, f"aucell_seed{seed}.pkl")
    with open(auc_file, 'wb') as f:
        pickle.dump(auc_matrix, f)
    print(f"Saved AUCell matrix: {auc_file}")
    
    print(f"pySCENIC workflow completed for seed {seed}")
    
    return adjacencies, modules, regulons, auc_matrix
