"""
GRNBoost Analysis Module for Organoid Morphogen Study

This module implements GRNBoost2-based network inference to identify
regulatory relationships between morphogens and transcription factors
in human brain organoid development.

Author: Research Team
Date: 2025
"""

import pandas as pd
import numpy as np
import pickle
import os
import re
from typing import List, Dict, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MorphogenNetworkAnalyzer:
    """
    Analyzes morphogen-regulon networks using GRNBoost2.
    
    This class combines AUCell regulon activity matrices with morphogen
    treatment metadata to infer regulatory networks between morphogens
    and transcription factors.
    """
    
    def __init__(self, cell_line: str, occur_threshold: int = 0):
        """
        Initialize the morphogen network analyzer.
        
        Args:
            cell_line: Cell line identifier (H1, H9, WIBJ2, WTC)
            occur_threshold: Occurrence threshold for consensus regulons
        """
        self.cell_line = cell_line
        self.occur_threshold = occur_threshold
        self.morphogens = ['FGF8', 'SHH', 'CHIR', 'RA']
        self.aucell_matrix = None
        self.metadata = None
        self.morphogen_matrix = None
        self.network_adjacencies = None
        
    def load_aucell_matrix(self, aucell_file: str) -> pd.DataFrame:
        """
        Load AUCell activity matrix from pickle file.
        
        Args:
            aucell_file: Path to AUCell pickle file
            
        Returns:
            AUCell activity matrix
        """
        logger.info(f"Loading AUCell matrix from {aucell_file}")
        
        with open(aucell_file, 'rb') as f:
            self.aucell_matrix = pickle.load(f)
            
        logger.info(f"Loaded AUCell matrix: {self.aucell_matrix.shape}")
        return self.aucell_matrix
    
    def prepare_morphogen_metadata(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare morphogen metadata for network analysis.
        
        Args:
            metadata: Cell metadata DataFrame
            
        Returns:
            Processed morphogen metadata
        """
        logger.info("Preparing morphogen metadata")
        
        meta_processed = metadata.copy()
        morphogen_vars = self.morphogens.copy()
        
        # Convert time to quantitative values
        if 'Time' in meta_processed.columns:
            meta_processed['Time_q'] = np.nan
            meta_processed.loc[meta_processed['Time'] == 'late', 'Time_q'] = 1
            meta_processed.loc[meta_processed['Time'] == 'middle', 'Time_q'] = 0
            meta_processed.loc[meta_processed['Time'] == 'early', 'Time_q'] = -1
            morphogen_vars.append('Time_q')
        
        # Convert medium to one-hot encoding
        if 'medium' in meta_processed.columns:
            for medium in set(meta_processed['medium']):
                meta_processed[medium] = 0
                meta_processed.loc[meta_processed['medium'] == medium, medium] = 1
                morphogen_vars.append(medium)
        
        # Convert time categories to one-hot encoding
        if 'Time' in meta_processed.columns:
            for time_point in set(meta_processed['Time']):
                col_name = f'timing_{time_point}'
                meta_processed[col_name] = 0
                meta_processed.loc[meta_processed['Time'] == time_point, col_name] = 1
                morphogen_vars.append(col_name)
        
        # Select morphogen variables
        available_vars = [var for var in morphogen_vars if var in meta_processed.columns]
        self.metadata = meta_processed[available_vars].copy()
        
        # Normalize morphogen variables
        self.metadata = np.log1p(self.metadata)
        self.metadata = (self.metadata - self.metadata.min()) / (self.metadata.max() - self.metadata.min())
        
        logger.info(f"Prepared morphogen metadata: {self.metadata.shape}")
        logger.info(f"Variables: {list(self.metadata.columns)}")
        
        return self.metadata
    
    def create_combined_matrix(self) -> pd.DataFrame:
        """
        Combine AUCell matrix with morphogen metadata for GRNBoost analysis.
        
        Returns:
            Combined matrix for network inference
        """
        logger.info("Creating combined AUCell + morphogen matrix")
        
        if self.aucell_matrix is None or self.metadata is None:
            raise ValueError("AUCell matrix and metadata must be loaded first")
        
        # Align cells between matrices
        common_cells = self.aucell_matrix.index.intersection(self.metadata.index)
        logger.info(f"Common cells: {len(common_cells)}")
        
        if len(common_cells) < 100:
            raise ValueError(f"Insufficient cell overlap: {len(common_cells)}")
        
        # Create combined matrix
        aucell_aligned = self.aucell_matrix.loc[common_cells]
        meta_aligned = self.metadata.loc[common_cells]
        
        # Add morphogens as "fake TFs" to the AUCell matrix
        for morph_var in meta_aligned.columns:
            aucell_aligned[morph_var] = meta_aligned[morph_var]
        
        self.morphogen_matrix = aucell_aligned
        
        logger.info(f"Combined matrix shape: {self.morphogen_matrix.shape}")
        logger.info(f"Morphogen variables added: {list(meta_aligned.columns)}")
        
        return self.morphogen_matrix
    
    def run_grnboost2(self, n_workers: int = 10, seed: int = 42) -> pd.DataFrame:
        """
        Run GRNBoost2 network inference on combined matrix.
        
        Args:
            n_workers: Number of workers for parallel processing
            seed: Random seed for reproducibility
            
        Returns:
            Network adjacency matrix
        """
        logger.info("Running GRNBoost2 network inference")
        
        try:
            from arboreto.algo import grnboost2
            from distributed import LocalCluster, Client
        except ImportError:
            raise ImportError("arboreto and dask packages required for GRNBoost2")
        
        if self.morphogen_matrix is None:
            raise ValueError("Combined matrix must be created first")
        
        # Define TF names (morphogen variables)
        tf_names = list(self.metadata.columns)
        logger.info(f"TF names (morphogens): {tf_names}")
        
        # Set up Dask cluster
        local_cluster = LocalCluster(
            n_workers=n_workers,
            threads_per_worker=2,
            memory_limit='8GB'
        )
        client = Client(local_cluster)
        
        try:
            # Run GRNBoost2
            logger.info("Starting GRNBoost2 inference...")
            self.network_adjacencies = grnboost2(
                expression_data=self.morphogen_matrix,
                tf_names=tf_names,
                client_or_address=client,
                seed=seed,
                verbose=True
            )
            
            logger.info(f"Generated {len(self.network_adjacencies)} network edges")
            
        finally:
            # Clean up Dask cluster
            client.close()
            local_cluster.close()
        
        return self.network_adjacencies
    
    def process_network_results(self, min_importance: float = 0.0) -> pd.DataFrame:
        """
        Process and filter network inference results.
        
        Args:
            min_importance: Minimum edge importance threshold
            
        Returns:
            Processed network summary
        """
        logger.info("Processing network results")
        
        if self.network_adjacencies is None:
            raise ValueError("Network adjacencies must be computed first")
        
        # Filter by importance threshold
        filtered_network = self.network_adjacencies[
            self.network_adjacencies['importance'] >= min_importance
        ].copy()
        
        # Add analysis metadata
        filtered_network['cell_line'] = self.cell_line
        filtered_network['occur_threshold'] = self.occur_threshold
        
        # Identify morphogen-regulon interactions
        morphogen_vars = list(self.metadata.columns)
        filtered_network['is_morphogen_tf'] = filtered_network['TF'].isin(morphogen_vars)
        filtered_network['is_morphogen_target'] = filtered_network['target'].isin(morphogen_vars)
        
        # Sort by importance
        filtered_network = filtered_network.sort_values('importance', ascending=False)
        
        logger.info(f"Filtered network: {len(filtered_network)} edges")
        logger.info(f"Morphogen TF edges: {filtered_network['is_morphogen_tf'].sum()}")
        
        return filtered_network
    
    def save_results(self, output_dir: str, save_adjacencies: bool = True,
                    save_modules: bool = True) -> None:
        """
        Save network analysis results.
        
        Args:
            output_dir: Output directory path
            save_adjacencies: Whether to save adjacency matrices
            save_modules: Whether to save module information
        """
        logger.info(f"Saving results to {output_dir}")
        
        os.makedirs(output_dir, exist_ok=True)
        
        if save_adjacencies and self.network_adjacencies is not None:
            # Save adjacencies
            adj_file = os.path.join(
                output_dir, 
                f'adjacencies_winteract_{self.occur_threshold}_combined.p'
            )
            with open(adj_file, 'wb') as f:
                pickle.dump(self.network_adjacencies, f)
            
            # Save as CSV for easy inspection
            csv_file = os.path.join(
                output_dir,
                f'adjacencies_winteract_{self.occur_threshold}_combined.csv'
            )
            self.network_adjacencies.to_csv(csv_file, index=False)
        
        if save_modules and self.morphogen_matrix is not None:
            # Create module summary similar to original analysis
            self._save_module_summary(output_dir)
    
    def _save_module_summary(self, output_dir: str) -> None:
        """
        Save module summary information.
        
        Args:
            output_dir: Output directory path
        """
        if self.network_adjacencies is None:
            return
        
        # Create module summary from network adjacencies
        module_summary = []
        
        for _, edge in self.network_adjacencies.iterrows():
            tf = edge['TF']
            target = edge['target']
            importance = edge['importance']
            
            # Skip self-interactions
            if tf == target:
                continue
            
            module_summary.append({
                'morph': tf,
                'gene': target,
                'w': importance,
                'context': 'morphogen_network'
            })
        
        module_df = pd.DataFrame(module_summary)
        
        # Save full summary
        summary_file = os.path.join(
            output_dir,
            f'module_summary_winteract_{self.occur_threshold}_combined.tsv'
        )
        module_df.to_csv(summary_file, sep='\t', index=False)
        
        # Create unique summary (highest weight per TF-gene pair)
        if len(module_df) > 0:
            module_df['morph-TF'] = module_df['morph'] + '_' + module_df['gene']
            module_unique = (module_df
                           .sort_values('w', ascending=False)
                           .groupby('morph-TF')
                           .head(1)
                           .copy())
            
            unique_file = os.path.join(
                output_dir,
                f'module_summary_uniq_{self.occur_threshold}_combined.tsv'
            )
            module_unique.to_csv(unique_file, sep='\t', index=False)
        
        logger.info(f"Saved module summaries: {len(module_df)} total edges")


def run_morphogen_network_analysis(cell_line: str, aucell_file: str, 
                                 metadata_file: str, output_dir: str,
                                 occur_threshold: int = 0,
                                 n_workers: int = 10) -> MorphogenNetworkAnalyzer:
    """
    Complete morphogen network analysis workflow.
    
    Args:
        cell_line: Cell line identifier
        aucell_file: Path to AUCell matrix file
        metadata_file: Path to metadata file
        output_dir: Output directory
        occur_threshold: Occurrence threshold for consensus regulons
        n_workers: Number of workers for GRNBoost2
        
    Returns:
        MorphogenNetworkAnalyzer instance with results
    """
    logger.info(f"Starting morphogen network analysis for {cell_line}")
    
    # Initialize analyzer
    analyzer = MorphogenNetworkAnalyzer(cell_line, occur_threshold)
    
    # Load data
    analyzer.load_aucell_matrix(aucell_file)
    
    # Load metadata (assuming it's an H5AD file)
    import scanpy as sc
    adata = sc.read_h5ad(metadata_file)
    analyzer.prepare_morphogen_metadata(adata.obs)
    
    # Create combined matrix
    analyzer.create_combined_matrix()
    
    # Run network inference
    analyzer.run_grnboost2(n_workers=n_workers)
    
    # Process results
    analyzer.process_network_results()
    
    # Save results
    analyzer.save_results(output_dir)
    
    logger.info(f"Morphogen network analysis completed for {cell_line}")
    
    return analyzer
