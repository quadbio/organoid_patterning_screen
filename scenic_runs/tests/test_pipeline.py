#!/usr/bin/env python3
"""
Test script for Organoid pySCENIC Pipeline
Tests the pipeline using existing pySCENIC results rather than running new analysis.
"""

import os
import sys
import logging
from pathlib import Path
import pickle
import pandas as pd
import scanpy as sc
import numpy as np

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from pyscenic_utils import load_config, load_existing_regulons, build_consensus_regulons
from morphogen_analysis import analyze_morphogen_correlations
from visualization import plot_regulon_summary, plot_morphogen_correlations

def setup_logging(config):
    """Setup logging configuration."""
    log_level = getattr(logging, config['logging']['level'].upper())
    log_file = config['logging']['log_file']
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    return logging.getLogger(__name__)

def test_load_existing_results(config, logger):
    """Test loading existing pySCENIC results."""
    logger.info("Testing loading of existing pySCENIC results...")
    
    # Get paths to existing results
    base_dir = Path(config['pyscenic']['existing_results']['base_dir'])
    seeds = config['pyscenic']['seeds']
    
    logger.info(f"Looking for results in: {base_dir}")
    logger.info(f"Testing with seeds: {seeds}")
    
    # Test loading regulons for each seed
    regulon_files = []
    for seed in seeds:
        regulon_file = base_dir / f"regulons_seed_{seed}.p"
        if regulon_file.exists():
            regulon_files.append(str(regulon_file))
            logger.info(f"Found regulon file: {regulon_file}")
        else:
            logger.warning(f"Missing regulon file: {regulon_file}")
    
    if not regulon_files:
        logger.error("No regulon files found!")
        return None
    
    # Load first regulon file to test
    try:
        with open(regulon_files[0], 'rb') as f:
            regulons = pickle.load(f)
        logger.info(f"Successfully loaded {len(regulons)} regulons from {regulon_files[0]}")
        
        # Show some regulon info
        for i, reg in enumerate(regulons[:3]):  # Show first 3 regulons
            logger.info(f"Regulon {i+1}: {reg.name} with {len(reg.gene2weight)} target genes")
            
        return regulon_files
        
    except Exception as e:
        logger.error(f"Error loading regulon file: {e}")
        return None

def test_consensus_building(config, regulon_files, logger):
    """Test building consensus regulons."""
    logger.info("Testing consensus regulon building...")
    
    try:
        # Load regulons from all seed files
        all_regulons = []
        for file_path in regulon_files[:5]:  # Test with first 5 seeds only
            with open(file_path, 'rb') as f:
                regulons = pickle.load(f)
                all_regulons.extend(regulons)
        
        logger.info(f"Loaded {len(all_regulons)} total regulons from {len(regulon_files[:5])} seed files")
        
        # Build consensus
        consensus_regulons = build_consensus_regulons(
            all_regulons, 
            min_seed_fraction=config['pyscenic']['consensus']['min_seed_fraction'],
            min_targets=config['pyscenic']['consensus']['min_targets']
        )
        
        logger.info(f"Built {len(consensus_regulons)} consensus regulons")
        
        # Show consensus regulon info
        for i, reg in enumerate(consensus_regulons[:3]):  # Show first 3
            logger.info(f"Consensus regulon {i+1}: {reg.name} with {len(reg.gene2weight)} target genes")
        
        return consensus_regulons
        
    except Exception as e:
        logger.error(f"Error building consensus regulons: {e}")
        return None

def test_data_loading(config, logger):
    """Test loading expression data."""
    logger.info("Testing expression data loading...")
    
    try:
        data_dir = Path(config['data']['input_dir'])
        count_files = config['data']['count_files']
        
        # Test loading one file
        test_file = data_dir / count_files['combined']
        logger.info(f"Loading test data from: {test_file}")
        
        if not test_file.exists():
            logger.error(f"Data file not found: {test_file}")
            return None
            
        adata = sc.read_h5ad(test_file)
        logger.info(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")
        
        return adata
        
    except Exception as e:
        logger.error(f"Error loading data: {e}")
        return None

def test_morphogen_analysis(config, adata, consensus_regulons, logger):
    """Test morphogen correlation analysis."""
    logger.info("Testing morphogen analysis...")
    
    try:
        morphogen_genes = config['analysis']['morphogens']['genes']
        logger.info(f"Analyzing morphogens: {morphogen_genes}")
        
        # Check which morphogen genes are present in data
        available_morphogens = [gene for gene in morphogen_genes if gene in adata.var_names]
        logger.info(f"Available morphogen genes in dataset: {available_morphogens}")
        
        if not available_morphogens:
            logger.warning("No morphogen genes found in dataset!")
            return None
            
        # For testing, create a simple correlation analysis
        correlations = analyze_morphogen_correlations(
            adata, 
            consensus_regulons, 
            available_morphogens,
            method=config['analysis']['correlation']['method']
        )
        
        logger.info(f"Computed correlations for {len(correlations)} regulon-morphogen pairs")
        
        return correlations
        
    except Exception as e:
        logger.error(f"Error in morphogen analysis: {e}")
        return None

def main():
    """Main test function."""
    # Load configuration
    config_file = "config/test_config.yaml"
    config = load_config(config_file)
    
    # Setup logging
    logger = setup_logging(config)
    logger.info("Starting pySCENIC pipeline test...")
    
    # Create output directory
    output_dir = Path(config['output']['results_dir'])
    output_dir.mkdir(exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Test 1: Load existing pySCENIC results
    regulon_files = test_load_existing_results(config, logger)
    if not regulon_files:
        logger.error("Failed to load existing results. Exiting.")
        return 1
    
    # Test 2: Build consensus regulons
    consensus_regulons = test_consensus_building(config, regulon_files, logger)
    if not consensus_regulons:
        logger.error("Failed to build consensus regulons. Exiting.")
        return 1
    
    # Test 3: Load expression data
    adata = test_data_loading(config, logger)
    if adata is None:
        logger.error("Failed to load expression data. Exiting.")
        return 1
    
    # Test 4: Morphogen analysis
    correlations = test_morphogen_analysis(config, adata, consensus_regulons, logger)
    
    # Save test results
    try:
        # Save consensus regulons
        consensus_file = output_dir / "consensus_regulons_test.p"
        with open(consensus_file, 'wb') as f:
            pickle.dump(consensus_regulons, f)
        logger.info(f"Saved consensus regulons to: {consensus_file}")
        
        # Save correlations if available
        if correlations is not None:
            corr_file = output_dir / "morphogen_correlations_test.csv"
            correlations.to_csv(corr_file)
            logger.info(f"Saved correlations to: {corr_file}")
        
        logger.info("Test completed successfully!")
        return 0
        
    except Exception as e:
        logger.error(f"Error saving results: {e}")
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
