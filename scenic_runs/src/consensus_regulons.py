"""
Consensus Regulon Generation for Organoid Morphogen Analysis

This module implements consensus regulon generation from multiple pySCENIC runs,
specifically designed for the organoid morphogen response analysis.

Author: Research Team
Date: 2025
"""

import os
import pickle
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Optional, Union
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ConsensusRegulonGenerator:
    """
    Generates consensus regulons from multiple pySCENIC runs for organoid analysis.
    
    This class processes multiple pySCENIC output files to create robust
    consensus regulons based on occurrence and size thresholds, following
    the methodology used in the organoid morphogen response study.
    """
    
    def __init__(self, sample_name: str, occur_threshold: int = 5, 
                 size_threshold: int = 10):
        """
        Initialize consensus regulon generator.
        
        Args:
            sample_name: Sample identifier (e.g., cell line name)
            occur_threshold: Minimum occurrences for TF-gene links
            size_threshold: Minimum regulon size threshold
        """
        self.sample_name = sample_name
        self.occur_threshold = occur_threshold
        self.size_threshold = size_threshold
        self.module_summary_df = None
        self.tf_gene_summary = None
        self.consensus_regulons = None
        
    def load_regulon_files(self, file_paths: List[str]) -> pd.DataFrame:
        """
        Load and summarize data from multiple regulon files.
        
        Args:
            file_paths: List of paths to regulon pickle files
            
        Returns:
            DataFrame with TF-gene relationships across runs
        """
        logger.info(f"Loading {len(file_paths)} regulon files for {self.sample_name}")
        
        module_summary = []
        
        for run_index, file_path in enumerate(file_paths, 1):
            logger.info(f"Processing run {run_index}: {file_path}")
            
            try:
                with open(file_path, 'rb') as file:
                    regulon_data = pickle.load(file)
                
                # Extract TF, gene, weight, and context from each regulon
                for regulon in regulon_data:
                    transcription_factor = re.sub("Regulon for ", "", regulon.name)
                    
                    for gene, weight in zip(regulon.genes, regulon.weights):
                        module_summary.append({
                            'TF': transcription_factor,
                            'gene': gene,
                            'weight': weight,
                            'context': regulon.context,
                            'run': str(run_index)
                        })
                        
            except Exception as e:
                logger.warning(f"Error processing {file_path}: {e}")
                continue
        
        self.module_summary_df = pd.DataFrame(module_summary)
        logger.info(f"Loaded data: {len(self.module_summary_df)} TF-gene relationships")
        
        return self.module_summary_df
    
    def plot_regulon_distribution(self) -> None:
        """
        Plot the distribution of detected regulons per run.
        """
        if self.module_summary_df is None:
            raise ValueError("Must load regulon files first")
        
        # Count regulons per run
        unique_runs = self.module_summary_df['run'].unique()
        regulon_counts = [
            len(set(self.module_summary_df[self.module_summary_df['run'] == run]['TF']))
            for run in unique_runs
        ]
        
        plt.figure(figsize=(8, 6))
        plt.boxplot(regulon_counts)
        plt.ylabel('Number of regulons detected', fontsize=12)
        plt.title(f'Regulon Detection Across Runs - {self.sample_name}')
        plt.show()
        
        logger.info(f"Regulon counts per run: {regulon_counts}")
    
    def calculate_tf_gene_statistics(self) -> pd.DataFrame:
        """
        Calculate occurrence and weight statistics for TF-gene links.
        
        Returns:
            DataFrame with TF-gene link statistics
        """
        if self.module_summary_df is None:
            raise ValueError("Must load regulon files first")
        
        logger.info("Calculating TF-gene link statistics")
        
        # Create unique TF-gene identifiers
        self.module_summary_df['TF_gene_link'] = (
            self.module_summary_df['TF'] + "__" + self.module_summary_df['gene']
        )
        
        # Initialize summary DataFrame
        self.tf_gene_summary = pd.DataFrame(
            index=self.module_summary_df['TF_gene_link'].unique()
        )
        
        # Calculate statistics
        self.tf_gene_summary['occurrences'] = (
            self.module_summary_df.groupby('TF_gene_link').size()
        )
        self.tf_gene_summary['mean_weight'] = (
            self.module_summary_df.groupby('TF_gene_link')['weight'].mean()
        )
        
        # Add TF and gene information
        first_occurrences = self.module_summary_df.groupby('TF_gene_link').first()
        self.tf_gene_summary['TF'] = first_occurrences['TF']
        self.tf_gene_summary['gene'] = first_occurrences['gene']
        
        logger.info(f"Calculated statistics for {len(self.tf_gene_summary)} TF-gene links")
        
        return self.tf_gene_summary
    
    def plot_occurrence_weight_distribution(self) -> None:
        """
        Plot the distribution of occurrences vs mean weights.
        """
        if self.tf_gene_summary is None:
            raise ValueError("Must calculate TF-gene statistics first")
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Violin plot
        sns.violinplot(x='occurrences', y='mean_weight', 
                      data=self.tf_gene_summary, ax=axes[0])
        axes[0].set_title('Occurrence vs Mean Weight Distribution')
        
        # Scatter plot
        axes[1].scatter(self.tf_gene_summary['occurrences'], 
                       self.tf_gene_summary['mean_weight'], 
                       c='black', alpha=0.6)
        axes[1].set_xlabel('Occurrences')
        axes[1].set_ylabel('Mean Weight')
        axes[1].set_title('Occurrence vs Mean Weight Scatter')
        
        plt.tight_layout()
        plt.show()
    
    def filter_high_quality_links(self) -> pd.DataFrame:
        """
        Filter TF-gene links based on occurrence and size thresholds.
        
        Returns:
            DataFrame with high-quality TF-gene links
        """
        if self.tf_gene_summary is None:
            raise ValueError("Must calculate TF-gene statistics first")
        
        logger.info("Filtering TF-gene links based on thresholds")
        logger.info(f"Occurrence threshold: {self.occur_threshold}")
        logger.info(f"Size threshold: {self.size_threshold}")
        
        # Filter by occurrence threshold
        high_quality_links = self.tf_gene_summary[
            self.tf_gene_summary['occurrences'] >= self.occur_threshold
        ].copy()
        
        # Filter by regulon size
        regulon_sizes = high_quality_links['TF'].value_counts()
        selected_regulons = regulon_sizes[regulon_sizes > self.size_threshold].index
        high_quality_links = high_quality_links[
            high_quality_links['TF'].isin(selected_regulons)
        ]
        
        logger.info(f"High-quality links: {len(high_quality_links)}")
        logger.info(f"Selected regulons: {len(selected_regulons)}")
        
        return high_quality_links
    
    def generate_consensus_regulons(self) -> List:
        """
        Generate consensus regulons from filtered TF-gene links.
        
        Returns:
            List of consensus regulon objects
        """
        try:
            import pyscenic as ps
        except ImportError:
            logger.warning("pySCENIC not available, creating mock regulons")
            return self._create_mock_regulons()
        
        high_quality_links = self.filter_high_quality_links()
        selected_tfs = high_quality_links['TF'].unique()
        
        logger.info("Creating consensus regulons")
        
        self.consensus_regulons = []
        
        for tf in selected_tfs:
            regulon_links = high_quality_links[high_quality_links['TF'] == tf]
            gene_weight_dict = dict(zip(regulon_links['gene'], regulon_links['mean_weight']))
            
            regulon = ps.utils.Regulon(
                name=tf,
                gene2weight=gene_weight_dict,
                gene2occurrence={},
                transcription_factor=re.sub(r"\(\+\)", "", tf)
            )
            self.consensus_regulons.append(regulon)
        
        logger.info(f"Generated {len(self.consensus_regulons)} consensus regulons")
        
        return self.consensus_regulons
    
    def _create_mock_regulons(self) -> List:
        """
        Create mock regulon objects when pySCENIC is not available.
        
        Returns:
            List of mock regulon objects
        """
        high_quality_links = self.filter_high_quality_links()
        selected_tfs = high_quality_links['TF'].unique()
        
        class MockRegulon:
            def __init__(self, name, genes, weights):
                self.name = name
                self.genes = genes
                self.weights = weights
                self.transcription_factor = re.sub(r"\(\+\)", "", name)
        
        self.consensus_regulons = []
        
        for tf in selected_tfs:
            regulon_links = high_quality_links[high_quality_links['TF'] == tf]
            genes = list(regulon_links['gene'])
            weights = list(regulon_links['mean_weight'])
            
            regulon = MockRegulon(tf, genes, weights)
            self.consensus_regulons.append(regulon)
        
        logger.info(f"Generated {len(self.consensus_regulons)} mock consensus regulons")
        
        return self.consensus_regulons
    
    def save_consensus_regulons(self, output_dir: str) -> None:
        """
        Save consensus regulons to files.
        
        Args:
            output_dir: Output directory path
        """
        if self.consensus_regulons is None:
            raise ValueError("Must generate consensus regulons first")
        
        logger.info(f"Saving consensus regulons to {output_dir}")
        
        # Create output directories
        sample_output_dir = os.path.join(output_dir, self.sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Save individual regulon gene lists
        for regulon in self.consensus_regulons:
            tf_cleaned_name = re.sub(r"\(\+\)", "", regulon.name)
            output_path = os.path.join(sample_output_dir, f"{tf_cleaned_name}.txt")
            
            with open(output_path, 'w') as file:
                file.write("\n".join(regulon.genes))
        
        # Save all consensus regulons as pickle
        consensus_pickle_path = os.path.join(output_dir, f"{self.sample_name}.p")
        with open(consensus_pickle_path, "wb") as file:
            pickle.dump(self.consensus_regulons, file)
        
        logger.info(f"Saved {len(self.consensus_regulons)} consensus regulons")
        logger.info(f"Individual gene lists: {sample_output_dir}/")
        logger.info(f"Combined pickle file: {consensus_pickle_path}")


def process_multiple_samples(sample_configs: Dict[str, Dict], 
                           output_base_dir: str = "regulons") -> Dict[str, List]:
    """
    Process consensus regulons for multiple samples in the organoid study.
    
    Args:
        sample_configs: Dictionary with sample configurations
        output_base_dir: Base output directory
        
    Returns:
        Dictionary mapping sample names to selected TFs
    """
    selected_tfs_all = {}
    
    for sample_name, config in sample_configs.items():
        logger.info(f"Processing sample: {sample_name}")
        
        generator = ConsensusRegulonGenerator(
            sample_name=sample_name,
            occur_threshold=config.get('occur_threshold', 5),
            size_threshold=config.get('size_threshold', 10)
        )
        
        # Load and process regulon files
        generator.load_regulon_files(config['file_paths'])
        generator.plot_regulon_distribution()
        generator.calculate_tf_gene_statistics()
        generator.plot_occurrence_weight_distribution()
        
        # Generate and save consensus regulons
        generator.generate_consensus_regulons()
        generator.save_consensus_regulons(output_base_dir)
        
        # Store selected TFs
        selected_tfs = [regulon.name for regulon in generator.consensus_regulons]
        selected_tfs_all[sample_name] = selected_tfs
        
        logger.info(f"Completed {sample_name}: {len(selected_tfs)} consensus regulons")
    
    return selected_tfs_all


def load_and_convert_motifs_to_regulons(motif_files: List[str], 
                                      output_dir: str) -> None:
    """
    Convert pySCENIC motif files to regulon pickle files.
    
    This function replicates the motif-to-regulon conversion process
    used in the original organoid morphogen analysis.
    
    Args:
        motif_files: List of enriched motifs CSV files
        output_dir: Output directory for regulon files
    """
    try:
        from pyscenic.utils import load_motifs
        from pyscenic.prune import df2regulons
    except ImportError:
        logger.error("pySCENIC not available for motif conversion")
        return
    
    logger.info(f"Converting {len(motif_files)} motif files to regulons")
    
    os.makedirs(output_dir, exist_ok=True)
    
    for file_path in motif_files:
        try:
            # Generate output filename
            new_file = re.sub("enriched_motifs_seed", "regulons_seed", file_path)
            new_file = re.sub(".csv", ".p", new_file)
            new_file = os.path.join(output_dir, os.path.basename(new_file))
            
            # Load motifs and convert to regulons
            df_motifs = load_motifs(file_path)
            regulons = df2regulons(df_motifs)
            
            # Save regulons
            with open(new_file, "wb") as f:
                pickle.dump(regulons, f)
                
            logger.info(f"Converted {file_path} -> {new_file}")
            
        except Exception as e:
            logger.warning(f"Error converting {file_path}: {e}")
            continue
    
    logger.info("Motif to regulon conversion completed")
