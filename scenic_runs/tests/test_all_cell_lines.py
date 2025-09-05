#!/usr/bin/env python3
"""
Comprehensive test of the organoid pySCENIC pipeline with multiple cell lines.
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

def test_all_cell_lines():
    """
    Test the pipeline with all available cell lines.
    """
    cell_lines = ['H1', 'H9', 'WIBJ2', 'WTC']
    morphogens = ['FGF8', 'SHH', 'CHIR', 'RA']
    results_summary = {}
    
    print("🧬 Testing Organoid pySCENIC Pipeline - Multi Cell Line Analysis")
    print("=" * 70)
    
    for cell_line in cell_lines:
        print(f"\n📊 Processing {cell_line}...")
        print("-" * 30)
        
        try:
            # Check for AUCell results
            aucell_file = f"../pyscenic_post/regulons/consensus_0/aucell_{cell_line}_combined.p"
            data_file = f"../data/exp1_counts_for_scenic_{cell_line}.h5ad"
            
            if not os.path.exists(aucell_file):
                print(f"❌ No AUCell file found for {cell_line}")
                continue
            if not os.path.exists(data_file):
                print(f"❌ No data file found for {cell_line}")
                continue
                
            # Load data
            with open(aucell_file, 'rb') as f:
                aucell_matrix = pickle.load(f)
            
            import scanpy as sc
            adata = sc.read_h5ad(data_file)
            
            print(f"✅ Data loaded - AUCell: {aucell_matrix.shape}, Cells: {adata.shape}")
            
            # Prepare metadata
            metadata = adata.obs.copy()
            morphogens_expanded = morphogens.copy()
            
            # Add time as numeric
            if 'Time' in metadata.columns:
                time_map = {'early': -1, 'middle': 0, 'late': 1}
                metadata['Time_numeric'] = metadata['Time'].map(time_map)
                morphogens_expanded.append('Time_numeric')
            
            # Add medium variables
            if 'medium' in metadata.columns:
                for medium in metadata['medium'].unique():
                    metadata[f'medium_{medium}'] = (metadata['medium'] == medium).astype(int)
                    morphogens_expanded.append(f'medium_{medium}')
            
            # Select available morphogens
            available_morphogens = [m for m in morphogens_expanded if m in metadata.columns]
            
            if not available_morphogens:
                print(f"❌ No morphogens found for {cell_line}")
                continue
            
            # Align data
            common_cells = aucell_matrix.index.intersection(metadata.index)
            if len(common_cells) < 100:
                print(f"❌ Insufficient overlap for {cell_line}: {len(common_cells)} cells")
                continue
            
            aucell_aligned = aucell_matrix.loc[common_cells]
            meta_aligned = metadata.loc[common_cells, available_morphogens]
            
            # Calculate correlations
            from scipy.stats import pearsonr
            
            correlations = []
            p_values = []
            
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
            
            corr_df = pd.DataFrame(correlations, 
                                  index=aucell_aligned.columns, 
                                  columns=available_morphogens)
            
            pval_df = pd.DataFrame(p_values, 
                                  index=aucell_aligned.columns, 
                                  columns=available_morphogens)
            
            # Store results
            results_summary[cell_line] = {
                'n_cells': len(common_cells),
                'n_regulons': len(aucell_aligned.columns),
                'n_morphogens': len(available_morphogens),
                'mean_abs_corr': corr_df.abs().mean().mean(),
                'max_corr': corr_df.max().max(),
                'min_corr': corr_df.min().min(),
                'n_significant': (pval_df < 0.05).sum().sum(),
                'morphogens': available_morphogens
            }
            
            # Save individual results
            output_dir = f'test_results_{cell_line}'
            os.makedirs(output_dir, exist_ok=True)
            corr_df.to_csv(f"{output_dir}/correlations.csv")
            pval_df.to_csv(f"{output_dir}/p_values.csv")
            
            print(f"✅ {cell_line} analysis complete:")
            print(f"   - Cells: {len(common_cells):,}")
            print(f"   - Regulons: {len(aucell_aligned.columns)}")
            print(f"   - Morphogens: {len(available_morphogens)}")
            print(f"   - Mean |correlation|: {corr_df.abs().mean().mean():.3f}")
            print(f"   - Significant correlations: {(pval_df < 0.05).sum().sum()}")
            
        except Exception as e:
            print(f"❌ Error processing {cell_line}: {e}")
            continue
    
    # Create summary report
    print("\n📋 SUMMARY REPORT")
    print("=" * 50)
    
    if results_summary:
        summary_df = pd.DataFrame(results_summary).T
        print(summary_df)
        
        # Save summary
        summary_df.to_csv('test_results/pipeline_summary.csv')
        
        # Create comparison plot
        plt.figure(figsize=(15, 10))
        
        # Plot 1: Number of significant correlations
        plt.subplot(2, 3, 1)
        plt.bar(summary_df.index, summary_df['n_significant'])
        plt.title('Significant Correlations by Cell Line')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        
        # Plot 2: Mean absolute correlation
        plt.subplot(2, 3, 2)
        plt.bar(summary_df.index, summary_df['mean_abs_corr'])
        plt.title('Mean Absolute Correlation')
        plt.ylabel('Correlation')
        plt.xticks(rotation=45)
        
        # Plot 3: Max correlation range
        plt.subplot(2, 3, 3)
        plt.bar(summary_df.index, summary_df['max_corr'], alpha=0.7, label='Max')
        plt.bar(summary_df.index, summary_df['min_corr'], alpha=0.7, label='Min')
        plt.title('Correlation Range')
        plt.ylabel('Correlation')
        plt.legend()
        plt.xticks(rotation=45)
        
        # Plot 4: Data dimensions
        plt.subplot(2, 3, 4)
        plt.bar(summary_df.index, summary_df['n_cells'])
        plt.title('Number of Cells')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        
        plt.subplot(2, 3, 5)
        plt.bar(summary_df.index, summary_df['n_regulons'])
        plt.title('Number of Regulons')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        
        plt.subplot(2, 3, 6)
        plt.bar(summary_df.index, summary_df['n_morphogens'])
        plt.title('Number of Morphogens')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        
        plt.tight_layout()
        plt.savefig('test_results/pipeline_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"\n✅ Pipeline test completed successfully!")
        print(f"Results saved to test_results/ directories")
        print(f"Summary saved to: test_results/pipeline_summary.csv")
        print(f"Comparison plot: test_results/pipeline_comparison.png")
        
    else:
        print("❌ No successful analyses completed")

if __name__ == "__main__":
    test_all_cell_lines()
