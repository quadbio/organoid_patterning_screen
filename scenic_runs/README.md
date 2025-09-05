# Organoid pySCENIC Pipeline

A complete 4-stage pipeline for pySCENIC regulatory network inference and morphogen analysis in organoids.

## 🔬 Overview

This pipeline implements a comprehensive workflow for analyzing regulatory networks in organoid development using pySCENIC. It includes SLURM-based parallel processing, consensus regulon generation, morphogen network inference, and correlation analysis.

## 📋 Table of Contents

- [Pipeline Stages](#pipeline-stages)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## 🔄 Pipeline Stages

### Stage 1: pySCENIC Runs
**Location**: `01_pyscenic_runs/`
**Purpose**: Run multiple pySCENIC analyses with different parameters
- Multi-array SLURM submission
- Stratified subsampling across conditions
- HVG selection and region-seed combinations
- Parallel processing across cell lines

### Stage 2: Consensus Regulons
**Location**: `02_consensus_regulons/`
**Purpose**: Generate consensus regulons from multiple runs
- Occurrence threshold filtering (20 for individual, 0 for combined)
- Size threshold filtering
- Quality control and validation

### Stage 3: Morphogen Networks
**Location**: `03_morphogen_networks/`
**Purpose**: Infer morphogen-regulon interaction networks
- GRNBoost2 network inference
- Morphogen → regulon connections
- Timing and medium condition analysis

### Stage 4: Final Analysis
**Location**: `04_final_analysis/`
**Purpose**: Calculate correlations and create publication figures
- Correlation matrix generation
- Statistical testing
- Publication-quality visualizations

## 🛠 Requirements

### System Requirements
- Linux/Unix system with SLURM scheduler
- Python 3.9+
- Conda/Mamba package manager
- Sufficient computational resources for parallel processing

### Python Dependencies
```
pyscenic==0.12.1
scanpy>=1.8.0
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
seaborn>=0.11.0
arboreto>=0.1.6
```

## 💾 Installation

```bash
# 1. Clone or create the pipeline directory
git clone <your-repo> organoid_pyscenic_pipeline
cd organoid_pyscenic_pipeline

# 2. Create conda environment
mamba create -n pyscenic python=3.9
mamba activate pyscenic

# 3. Install dependencies
pip install pyscenic==0.12.1
pip install scanpy pandas numpy matplotlib seaborn arboreto

# 4. Verify installation
python -c "import pyscenic; print(f'pySCENIC {pyscenic.__version__} installed')"
```

## 🚀 Quick Start

1. **Prepare the data** (AnnData format in `data/`):
   ```
   data/
   ├── exp1_counts_for_scenic_H1.h5ad
   ├── exp1_counts_for_scenic_WTC.h5ad
   ├── exp1_counts_for_scenic_H9.h5ad
   └── exp1_counts_for_scenic_WIBJ2.h5ad
   ```

2. **Run Stage 1** (pySCENIC runs):
   ```bash
   cd 01_pyscenic_runs
   
   # Generate parameter combinations
   jupyter notebook generate_combinations.ipynb
   
   # Submit SLURM array jobs
   sbatch submit_multi_pyscenic.sh
   ```

3. **Run Stage 2** (consensus regulons):
   ```bash
   cd ../02_consensus_regulons
   jupyter notebook consensus_generation.ipynb
   ```

4. **Run Stage 3** (morphogen networks):
   ```bash
   cd ../03_morphogen_networks
   jupyter notebook morphogen_analysis.ipynb
   ```

5. **Run Stage 4** (final analysis):
   ```bash
   cd ../04_final_analysis
   jupyter notebook correlation_analysis.ipynb
   ```

## 📖 Detailed Usage

### Stage 1: pySCENIC Runs

**Key Files**:
- `generate_combinations.ipynb`: Create region-seed parameter combinations
- `run_multi_pyscenic.py`: Main pySCENIC execution script
- `submit_multi_pyscenic.sh`: SLURM submission script

**Utilities**: Core functions are in `src/pyscenic_utils.py` for data processing

**Parameters**:
- **Subsampling**: 5000 cells per condition (stratified)
- **HVG selection**: 2000 most variable genes
- **Seeds**: Multiple random seeds for robustness
- **Regions**: Different TF binding site databases

**Execution**:
```bash
# Generate 100 region-seed combinations
python -c "from pyscenic_utils import generate_combinations; generate_combinations(n_combinations=100)"

# Submit array job (modify SLURM parameters as needed)
sbatch submit_multi_pyscenic.sh
```

### Stage 2: Consensus Regulons

**Key Parameters**:
- **Individual cell lines**: `occur_threshold=20`, `size_threshold=0`
- **Combined analysis**: `occur_threshold=0`, `size_threshold=0`

**Process**:
1. Load regulon files from multiple pySCENIC runs
2. Calculate occurrence frequency across runs
3. Filter by thresholds
4. Generate consensus regulon sets

### Stage 3: Morphogen Networks

**Method**: GRNBoost2 network inference
**Input**: Morphogen/timing data + AUCell regulon activities
**Output**: Morphogen → regulon interaction networks

**Key Features**:
- Connects morphogens to regulon activities
- Includes timing and medium conditions
- Generates importance scores for interactions

### Stage 4: Final Analysis

**Outputs**:
- Correlation matrices (CSV)
- Publication-quality heatmaps (PNG)
- Network visualizations (PNG)
- Summary statistics (CSV)

## 📁 Output Files

```
organoid_pyscenic_pipeline/
├── 01_pyscenic_runs/
│   ├── results/
│   │   ├── H1/
│   │   │   ├── regulons_seed*.pkl
│   │   │   └── aucell_seed*.csv
│   │   └── [WTC, H9, WIBJ2]/
│   └── region_seed_combos.txt
├── 02_consensus_regulons/
│   └── regulons/
│       ├── consensus_20/  # Individual cell lines
│       └── consensus_0/   # Combined analysis
├── 03_morphogen_networks/
│   └── networks/
│       ├── morphogen_regulon_network_*.csv
│       └── morphogen_regulon_networks_combined.csv
└── 04_final_analysis/
    ├── final_correlations_combined.csv
    ├── correlation_matrix.csv
    ├── summary_statistics.csv
    └── plots/
        ├── correlation_heatmap.png
        ├── correlation_distribution.png
        └── network_graph.png
```

## 🔧 Troubleshooting

### Common Issues

1. **pySCENIC version conflicts**:
   ```bash
   pip install pyscenic==0.12.1 --force-reinstall
   ```

2. **Memory issues with large datasets**:
   - Increase SLURM memory allocation
   - Reduce subsampling size in `pyscenic_utils.py`

3. **SLURM job failures**:
   - Check `pyscenic_logs/` for error messages
   - Verify resource allocations in submission script

4. **Missing morphogen data**:
   - Ensure morphogen metadata is available
   - Check `load_morphogen_data()` function in `grnboost_analysis.py`

### Performance Optimization

- **Parallel processing**: Adjust `--n_workers` in pySCENIC calls
- **Memory usage**: Monitor with `squeue` and adjust accordingly
- **Storage**: Use fast local storage for temporary files

## 📊 Expected Results

- **Stage 1**: ~100 regulon files per cell line
- **Stage 2**: Consensus regulons (typically 50-200 per cell line)
- **Stage 3**: Morphogen-regulon networks (hundreds of interactions)
- **Stage 4**: Publication-ready correlation matrices and visualizations

## 📝 Citation

Please cite the following in your publications:

**Primary reference**:
```
Decoding morphogen patterning of human neural organoids with a multiplexed single-cell transcriptomic screen
Fátima Sanchís-Calleja, Akanksha Jain, Zhisong He, Ryoko Okamoto, Charlotte Rusimbi, Pedro Rifes, Gaurav Singh Rathore, Malgorzata Santel, Jasper Janssens, Makiko Seimiya, Jonas Simon Fleck, Agnete Kirkeby, J. Gray Camp, Barbara Treutlein
https://www.biorxiv.org/content/10.1101/2024.02.08.579413v1.full
```

**Methods**:
```
pySCENIC: Aibar et al. (2017). SCENIC: single-cell regulatory network inference and clustering. Nature Methods.
GRNBoost2: Moerman et al. (2019). GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks.
```

---

**Last updated**: March 2025
**Version**: 1.0.0
**Compatible with**: pySCENIC 0.12.1
