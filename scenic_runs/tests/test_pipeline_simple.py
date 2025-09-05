#!/usr/bin/env python3
"""
Test script for the organoid pySCENIC pipeline using existing results.
This version avoids problematic pySCENIC imports and focuses on morphogen analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add src directory to path
sys.path.insert(0, 'src')

def load_existing_aucell_results(results_dir, cell_line='H1'):
    """
    Load existing AUCell results from the original analysis.
    """
    print(f"Loading AUCell results for {cell_line}...")
    
    # Try different possible locations for AUCell results
    possible_files = [
        f"../pyscenic_post/regulons/consensus_0/aucell_{cell_line}_combined.p",
        f"../pyscenic_celline/results/{cell_line}/aucell_matrix.p",
        f"results/aucell_{cell_line}_combined.p"
    ]
    
    for file_path in possible_files:
        if os.path.exists(file_path):
            print(f"Found AUCell file: {file_path}")
            with open(file_path, 'rb') as f:
                aucell_matrix = pickle.load(f)
            print(f"AUCell matrix shape: {aucell_matrix.shape}")
            return aucell_matrix
    
    # If no real file found, create a mock dataset for testing
    print("No AUCell file found, creating mock data for testing...")
    n_cells = 1000
    n_regulons = 50
    np.random.seed(42)
    
    # Create mock regulon names
    regulon_names = [f"TF{i}_regulon" for i in range(n_regulons)]
    cell_names = [f"cell_{i}" for i in range(n_cells)]
    
    # Create mock AUCell matrix with some structure
    aucell_matrix = pd.DataFrame(
        np.random.beta(2, 5, (n_cells, n_regulons)),  # Beta distribution for AUCell-like values
        index=cell_names,
        columns=regulon_names
    )
    
    return aucell_matrix

def load_sample_data(cell_line='H1'):
    """
    Load sample data for the specified cell line.
    """
    print(f"Loading sample data for {cell_line}...")
    
    # Try to find the actual data file
    data_files = [
        f"../data/exp1_counts_for_scenic_{cell_line}.h5ad",
        f"data/exp1_counts_for_scenic_{cell_line}.h5ad"
    ]
    
    for data_file in data_files:
        if os.path.exists(data_file):
            print(f"Found data file: {data_file}")
            try:
                import scanpy as sc
                adata = sc.read_h5ad(data_file)
                print(f"Data shape: {adata.shape}")
                return adata
            except Exception as e:
                print(f"Error loading data: {e}")
                
    # Create mock metadata for testing
    print("Creating mock metadata for testing...")
    n_cells = 1000
    np.random.seed(42)
    
    # Create mock metadata that matches the AUCell matrix
    metadata = pd.DataFrame({
        'cell_line': [cell_line] * n_cells,
        'FGF8': np.random.choice([0, 1, 2], n_cells),
        'SHH': np.random.choice([0, 1, 2], n_cells),
        'CHIR': np.random.choice([0, 1, 2], n_cells),
        'RA': np.random.choice([0, 1, 2], n_cells),
        'Time': np.random.choice(['early', 'middle', 'late'], n_cells),
        'medium': np.random.choice(['NPM', 'NIM'], n_cells),
        'RNA_snn_res.1': np.random.choice(range(20), n_cells)
    }, index=[f"cell_{i}" for i in range(n_cells)])
    
    # Create a mock AnnData-like object
    class MockAnnData:
        def __init__(self, metadata):
            self.obs = metadata
            self.shape = (len(metadata), 2000)  # Mock shape
    
    return MockAnnData(metadata)

def simple_morphogen_correlation_analysis(aucell_matrix, metadata, morphogens):
    """
    Simplified morphogen correlation analysis without pySCENIC dependencies.
    """
    print("Running morphogen correlation analysis...")
    
    # Prepare morphogen metadata
    meta_processed = metadata.copy()
    
    # Convert time to numeric
    if 'Time' in meta_processed.columns:
        time_map = {'early': -1, 'middle': 0, 'late': 1}
        meta_processed['Time_numeric'] = meta_processed['Time'].map(time_map)
        morphogens.append('Time_numeric')
    
    # One-hot encode categorical variables
    if 'medium' in meta_processed.columns:
        for medium in meta_processed['medium'].unique():
            meta_processed[f'medium_{medium}'] = (meta_processed['medium'] == medium).astype(int)
            morphogens.append(f'medium_{medium}')
    
    # Select morphogen columns that exist
    available_morphogens = [m for m in morphogens if m in meta_processed.columns]
    print(f"Available morphogens: {available_morphogens}")
    
    if not available_morphogens:
        print("No morphogens found in metadata!")
        return None, None
    
    # Align cells between AUCell matrix and metadata
    common_cells = aucell_matrix.index.intersection(meta_processed.index)
    print(f"Common cells: {len(common_cells)}")
    
    if len(common_cells) < 10:
        print("Not enough common cells for analysis!")
        return None, None
    
    aucell_aligned = aucell_matrix.loc[common_cells]
    meta_aligned = meta_processed.loc[common_cells, available_morphogens]
    
    # Calculate correlations
    correlations = []
    p_values = []
    
    from scipy.stats import pearsonr
    
    for regulon in aucell_aligned.columns:
        regulon_corrs = []
        regulon_pvals = []
        
        for morphogen in available_morphogens:
            try:
                corr, pval = pearsonr(aucell_aligned[regulon], meta_aligned[morphogen])
                regulon_corrs.append(corr)
                regulon_pvals.append(pval)
            except:
                regulon_corrs.append(0)
                regulon_pvals.append(1)
        
        correlations.append(regulon_corrs)
        p_values.append(regulon_pvals)
    
    # Create result DataFrames
    corr_df = pd.DataFrame(correlations, 
                          index=aucell_aligned.columns, 
                          columns=available_morphogens)
    
    pval_df = pd.DataFrame(p_values, 
                          index=aucell_aligned.columns, 
                          columns=available_morphogens)
    
    return corr_df, pval_df

def create_simple_plots(corr_df, pval_df, output_dir):
    """
    Create simple visualization plots.
    """
    print("Creating visualization plots...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Correlation heatmap
    plt.figure(figsize=(10, 8))
    
    # Select top regulons by max absolute correlation
    max_abs_corr = corr_df.abs().max(axis=1)
    top_regulons = max_abs_corr.nlargest(20).index
    
    plot_data = corr_df.loc[top_regulons]
    significance_mask = (pval_df.loc[top_regulons] < 0.05)
    
    sns.heatmap(plot_data, 
                annot=False,
                cmap='RdBu_r',
                center=0,
                vmin=-1, vmax=1,
                cbar_kws={'label': 'Correlation'})
    
    # Add significance markers
    for i, regulon in enumerate(top_regulons):
        for j, morphogen in enumerate(plot_data.columns):
            if significance_mask.loc[regulon, morphogen]:
                plt.text(j + 0.5, i + 0.5, '*', 
                        ha='center', va='center', 
                        color='white', fontsize=8, fontweight='bold')
    
    plt.title('Top 20 Regulon-Morphogen Correlations\n(* p < 0.05)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/correlation_heatmap.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Correlation distribution
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.hist(corr_df.values.flatten(), bins=50, alpha=0.7)
    plt.xlabel('Correlation')
    plt.ylabel('Frequency')
    plt.title('Correlation Distribution')
    
    plt.subplot(1, 3, 2)
    plt.hist(max_abs_corr.values, bins=30, alpha=0.7)
    plt.xlabel('Max Absolute Correlation')
    plt.ylabel('Frequency')
    plt.title('Max Correlation per Regulon')
    
    plt.subplot(1, 3, 3)
    significant_corrs = corr_df[pval_df < 0.05]
    plt.bar(range(len(corr_df.columns)), (pval_df < 0.05).sum(), alpha=0.7)
    plt.xticks(range(len(corr_df.columns)), corr_df.columns, rotation=45)
    plt.ylabel('Number of Significant Correlations')
    plt.title('Significant Correlations by Morphogen')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/correlation_summary.png", dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """
    Main test function.
    """
    print("🧬 Testing Organoid pySCENIC Pipeline")
    print("=" * 50)
    
    # Configuration
    cell_line = 'H1'
    morphogens = ['FGF8', 'SHH', 'CHIR', 'RA']
    output_dir = 'test_results'
    
    try:
        # 1. Load existing AUCell results
        aucell_matrix = load_existing_aucell_results('results', cell_line)
        
        # 2. Load sample data
        adata = load_sample_data(cell_line)
        
        # 3. Run morphogen correlation analysis
        corr_df, pval_df = simple_morphogen_correlation_analysis(
            aucell_matrix, adata.obs, morphogens.copy()
        )
        
        if corr_df is not None:
            # 4. Print summary statistics
            print("\nAnalysis Results:")
            print("-" * 30)
            print(f"Correlation matrix shape: {corr_df.shape}")
            print(f"Mean absolute correlation: {corr_df.abs().mean().mean():.3f}")
            print(f"Max correlation: {corr_df.max().max():.3f}")
            print(f"Min correlation: {corr_df.min().min():.3f}")
            print(f"Significant correlations (p<0.05): {(pval_df < 0.05).sum().sum()}")
            
            # 5. Create visualizations
            create_simple_plots(corr_df, pval_df, output_dir)
            
            # 6. Save results
            corr_df.to_csv(f"{output_dir}/correlations.csv")
            pval_df.to_csv(f"{output_dir}/p_values.csv")
            
            print(f"\n✅ Test completed successfully!")
            print(f"Results saved to: {output_dir}/")
            print(f"Plots created: correlation_heatmap.png, correlation_summary.png")
            
        else:
            print("❌ Analysis failed - no results generated")
            
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
