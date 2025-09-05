"""
Correlation Analysis Utilities

This module provides functions for calculating correlations between morphogens and regulons,
creating correlation matrices, and generating publication-quality visualizations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os


def calculate_correlations(network_df, min_importance=0.1):
    """
    Calculate correlations from network importance scores.
    
    This method uses network importance as a proxy for correlation strength,
    following the GRNBoost2 methodology where importance scores reflect
    the strength of regulatory relationships.
    
    Parameters
    ----------
    network_df : pd.DataFrame
        Network dataframe with columns 'TF', 'target', 'importance', 'cell_line'
    min_importance : float, default=0.1
        Minimum importance threshold for filtering significant interactions
    
    Returns
    -------
    pd.DataFrame
        Correlation dataframe with columns 'morphogen', 'regulon', 'correlation', 'cell_line'
    """
    
    # Filter by minimum importance
    significant = network_df[network_df['importance'] >= min_importance].copy()
    
    # Create correlation entries
    correlations = []
    
    for _, row in significant.iterrows():
        correlations.append({
            'morphogen': row['TF'],
            'regulon': row['target'],
            'correlation': row['importance'],  # Using importance as correlation
            'cell_line': row.get('cell_line', 'unknown')
        })
    
    return pd.DataFrame(correlations)


def create_correlation_matrix(correlations_df):
    """
    Create a correlation matrix from correlation data.
    
    Parameters
    ----------
    correlations_df : pd.DataFrame
        Correlation dataframe with columns 'morphogen', 'regulon', 'correlation'
    
    Returns
    -------
    pd.DataFrame
        Pivot table with morphogens as rows and regulons as columns
    """
    # Pivot to create matrix
    matrix = correlations_df.pivot_table(
        index='morphogen', 
        columns='regulon', 
        values='correlation',
        fill_value=0
    )
    
    return matrix


def plot_correlation_heatmap(corr_matrix, top_n=20, figsize=(12, 8), output_path=None):
    """
    Create a heatmap of the top correlations.
    
    Parameters
    ----------
    corr_matrix : pd.DataFrame
        Correlation matrix with morphogens as rows and regulons as columns
    top_n : int, default=20
        Number of top morphogens and regulons to show
    figsize : tuple, default=(12, 8)
        Figure size
    output_path : str, optional
        Path to save the plot
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Only show top correlations for clarity
    top_morphogens = corr_matrix.max(axis=1).nlargest(top_n).index
    top_regulons = corr_matrix.max(axis=0).nlargest(top_n).index
    
    subset_matrix = corr_matrix.loc[top_morphogens, top_regulons]
    
    sns.heatmap(subset_matrix, cmap='RdYlBu_r', center=0, 
                cbar_kws={'label': 'Correlation'}, 
                xticklabels=True, yticklabels=True, ax=ax)
    
    ax.set_title('Top Morphogen-Regulon Correlations')
    ax.set_xlabel('Regulons')
    ax.set_ylabel('Morphogens')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_correlation_distribution(correlations_df, figsize=(10, 6), output_path=None):
    """
    Plot distribution of correlations by cell line.
    
    Parameters
    ----------
    correlations_df : pd.DataFrame
        Correlation dataframe with 'cell_line' and 'correlation' columns
    figsize : tuple, default=(10, 6)
        Figure size
    output_path : str, optional
        Path to save the plot
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Box plot of correlations by cell line
    sns.boxplot(data=correlations_df, x='cell_line', y='correlation', ax=ax)
    ax.set_title('Distribution of Morphogen-Regulon Correlations by Cell Line')
    ax.set_xlabel('Cell Line')
    ax.set_ylabel('Correlation Strength')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_network_graph(correlations_df, top_n=50, top_edges=20, figsize=(12, 10), output_path=None):
    """
    Create a network graph visualization of top interactions.
    
    Parameters
    ----------
    correlations_df : pd.DataFrame
        Correlation dataframe
    top_n : int, default=50
        Number of top interactions to consider
    top_edges : int, default=20
        Number of top edges to display
    figsize : tuple, default=(12, 10)
        Figure size
    output_path : str, optional
        Path to save the plot
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get top interactions
    top_interactions = correlations_df.nlargest(top_n, 'correlation')
    
    # Get unique nodes
    unique_morphogens = top_interactions['morphogen'].unique()
    unique_regulons = top_interactions['regulon'].unique()
    
    # Plot morphogens on left, regulons on right
    morph_y = np.linspace(0, 1, len(unique_morphogens))
    reg_y = np.linspace(0, 1, len(unique_regulons))
    
    # Plot nodes
    for i, morph in enumerate(unique_morphogens):
        ax.scatter(0, morph_y[i], s=100, c='red', alpha=0.7)
        ax.text(-0.05, morph_y[i], morph, ha='right', va='center', fontsize=8)
    
    for i, reg in enumerate(unique_regulons):
        ax.scatter(1, reg_y[i], s=100, c='blue', alpha=0.7)
        ax.text(1.05, reg_y[i], reg, ha='left', va='center', fontsize=8)
    
    # Plot edges
    for _, row in top_interactions.head(top_edges).iterrows():
        morph_idx = list(unique_morphogens).index(row['morphogen'])
        reg_idx = list(unique_regulons).index(row['regulon'])
        
        ax.plot([0, 1], [morph_y[morph_idx], reg_y[reg_idx]], 
               'k-', alpha=row['correlation'], linewidth=2*row['correlation'])
    
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title('Top Morphogen-Regulon Network Interactions')
    ax.set_xlabel('Morphogens → Regulons')
    ax.axis('off')
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    return fig


def calculate_summary_statistics(correlations_dict):
    """
    Calculate summary statistics for correlation results.
    
    Parameters
    ----------
    correlations_dict : dict
        Dictionary with cell_line as keys and correlation DataFrames as values
    
    Returns
    -------
    pd.DataFrame
        Summary statistics with one row per cell line
    """
    summary_stats = []
    
    for cell_line, corr_df in correlations_dict.items():
        stats_dict = {
            'cell_line': cell_line,
            'total_interactions': len(corr_df),
            'unique_morphogens': corr_df['morphogen'].nunique(),
            'unique_regulons': corr_df['regulon'].nunique(),
            'mean_correlation': corr_df['correlation'].mean(),
            'max_correlation': corr_df['correlation'].max(),
            'min_correlation': corr_df['correlation'].min(),
            'std_correlation': corr_df['correlation'].std()
        }
        summary_stats.append(stats_dict)
    
    return pd.DataFrame(summary_stats)


def save_correlation_results(correlations_combined, correlations_individual, 
                           corr_matrix, summary_stats, output_dir="."):
    """
    Save all correlation results to files.
    
    Parameters
    ----------
    correlations_combined : pd.DataFrame
        Combined correlations across all cell lines
    correlations_individual : dict
        Individual cell line correlations
    corr_matrix : pd.DataFrame
        Correlation matrix
    summary_stats : pd.DataFrame
        Summary statistics
    output_dir : str, default="."
        Output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Save combined correlations
    if correlations_combined is not None:
        correlations_combined.to_csv(f'{output_dir}/final_correlations_combined.csv', index=False)
        print("✅ Saved combined correlations")
    
    # Save correlation matrix
    if corr_matrix is not None:
        corr_matrix.to_csv(f'{output_dir}/correlation_matrix.csv')
        print("✅ Saved correlation matrix")
    
    # Save individual cell line correlations
    for cell_line, corr_df in correlations_individual.items():
        corr_df.to_csv(f'{output_dir}/final_correlations_{cell_line}.csv', index=False)
        print(f"✅ Saved correlations for {cell_line}")
    
    # Save summary statistics
    summary_stats.to_csv(f'{output_dir}/summary_statistics.csv', index=False)
    print("✅ Saved summary statistics")
