#!/bin/bash
#SBATCH --job-name=pyscenic_multi
#SBATCH --output=logs/pyscenic_%A_%a.out
#SBATCH --error=logs/pyscenic_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=20
#SBATCH --array=1-1000%50

# pySCENIC Multi-Array Job Submission Script
# Based on original pyscenic_celline submission

# Create logs directory
mkdir -p logs

# Load required modules
module load gcc/8.2.0 python_gpu/3.9.9

# Activate environment
source ../../.venv/bin/activate

# Configuration file
CONFIG_FILE="../config/pyscenic_config.yaml"

# Read region-seed combinations
COMBOS_FILE="region_seed_combos.txt"
COMBO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $COMBOS_FILE)

if [ -z "$COMBO" ]; then
    echo "No combination found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Parse combination
REGION=$(echo $COMBO | cut -d' ' -f1)
SEED=$(echo $COMBO | cut -d' ' -f2)

echo "Processing combination $SLURM_ARRAY_TASK_ID: Region=$REGION, Seed=$SEED"

# Set output directory
OUTPUT_DIR="results"

# Run the analysis
if [ "$REGION" = "cell_line" ]; then
    # Run for each cell line individually
    for CELL_LINE in H1 H9 WIBJ2 WTC; do
        echo "Running pySCENIC for cell line: $CELL_LINE, seed: $SEED"
        python run_multi_pyscenic.py \
            --config $CONFIG_FILE \
            --region cell_line \
            --cell_line $CELL_LINE \
            --seed $SEED \
            --output_dir $OUTPUT_DIR
        
        if [ $? -ne 0 ]; then
            echo "Error: pySCENIC failed for $CELL_LINE with seed $SEED"
            exit 1
        fi
    done
else
    # Run combined analysis
    echo "Running combined pySCENIC analysis, seed: $SEED"
    python run_multi_pyscenic.py \
        --config $CONFIG_FILE \
        --region combined \
        --seed $SEED \
        --output_dir $OUTPUT_DIR
    
    if [ $? -ne 0 ]; then
        echo "Error: Combined pySCENIC failed with seed $SEED"
        exit 1
    fi
fi

echo "✅ Completed pySCENIC analysis for combination $SLURM_ARRAY_TASK_ID"
