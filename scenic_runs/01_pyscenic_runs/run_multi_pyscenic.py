#!/usr/bin/env python3
"""
Multi-pySCENIC Runner for Organoid Analysis
Adapted from pyscenic_celline/multi_pyscenic.py

This script runs multiple pySCENIC analyses with different seeds and parameters
for both individual cell lines and combined stratified analysis.
"""

import argparse
import os
import sys
import yaml
from pathlib import Path
import pandas as pd
import scanpy as sc

# Add utils to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from pyscenic_utils import run_pyscenic_workflow

def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_expression_data(data_path, cell_line=None):
    """Load expression data for specific cell line or combined."""
    if cell_line:
        file_path = f"{data_path}/exp1_counts_for_scenic_{cell_line}.h5ad"
    else:
        file_path = f"{data_path}/exp1_counts_for_scenic.h5ad"
    
    print(f"Loading data from: {file_path}")
    adata = sc.read_h5ad(file_path)
    print(f"Loaded data shape: {adata.shape}")
    
    return adata

def main():
    parser = argparse.ArgumentParser(description='Run multi-pySCENIC analysis')
    parser.add_argument('--config', type=str, required=True,
                       help='Configuration YAML file')
    parser.add_argument('--region', type=str, required=True,
                       choices=['cell_line', 'combined'],
                       help='Analysis region: individual cell line or combined')
    parser.add_argument('--cell_line', type=str, 
                       choices=['H1', 'H9', 'WIBJ2', 'WTC'],
                       help='Cell line (required if region=cell_line)')
    parser.add_argument('--seed', type=int, required=True,
                       help='Random seed for this run')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.region == 'cell_line' and not args.cell_line:
        parser.error("--cell_line is required when region=cell_line")
    
    # Load configuration
    config = load_config(args.config)
    
    # Set up paths
    data_path = config['data']['input_dir']
    tf_file = config['pyscenic']['tf_file']
    ranking_dbs = config['pyscenic']['ranking_dbs']
    motif_annotations = config['pyscenic']['motif_annotations']
    
    # Load TF names
    tf_names = []
    with open(tf_file, 'r') as f:
        tf_names = [line.strip() for line in f if line.strip()]
    
    print(f"Loaded {len(tf_names)} transcription factors")
    
    # Load expression data
    if args.region == 'cell_line':
        adata = load_expression_data(data_path, args.cell_line)
        output_subdir = os.path.join(args.output_dir, args.cell_line)
        subsample_column = config['subsampling'].get('column', 'cell_line')
    else:
        adata = load_expression_data(data_path)
        output_subdir = os.path.join(args.output_dir, 'combined')
        subsample_column = config['subsampling'].get('stratify_column', 'cell_line')
    
    os.makedirs(output_subdir, exist_ok=True)
    
    # Get pySCENIC parameters
    pyscenic_params = config['pyscenic']
    subsample_params = config['subsampling']
    
    # Run pySCENIC workflow
    try:
        adjacencies, modules, regulons, auc_matrix = run_pyscenic_workflow(
            adata=adata,
            tf_names=tf_names,
            ranking_dbs=ranking_dbs,
            motif_annotations=motif_annotations,
            output_dir=output_subdir,
            seed=args.seed,
            hvg_genes=pyscenic_params.get('hvg_genes', 3000),
            subsample_column=subsample_column,
            subsample_ncells=subsample_params.get('ncells'),
            strat=subsample_params.get('strategy'),
            n_workers=pyscenic_params.get('n_workers', 8),
            memory_limit=pyscenic_params.get('memory_limit', "8e9")
        )
        
        print(f"✅ pySCENIC completed successfully for {args.region}")
        if args.cell_line:
            print(f"   Cell line: {args.cell_line}")
        print(f"   Seed: {args.seed}")
        print(f"   Output: {output_subdir}")
        
    except Exception as e:
        print(f"❌ pySCENIC failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
