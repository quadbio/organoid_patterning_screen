# Data Directory

This directory should contain your input data files for the pySCENIC pipeline.

## Expected Files

Place your AnnData files here with the following naming convention:

```
data/
├── exp1_counts_for_scenic_H1.h5ad
├── exp1_counts_for_scenic_WTC.h5ad  
├── exp1_counts_for_scenic_H9.h5ad
└── exp1_counts_for_scenic_WIBJ2.h5ad
```

## Data Format Requirements

Each H5AD file should contain:
- **Raw counts** in `.X` 
- **Cell metadata** in `.obs` with columns for:
  - Cell line information
  - Morphogen treatment conditions  
  - Time points
  - Medium conditions
- **Gene names** in `.var_names`

## Data Preparation

If your data is in a different format, use scanpy to convert:

```python
import scanpy as sc
import pandas as pd

# Load your data
adata = sc.read_csv("your_counts.csv")  # or other format
meta = pd.read_csv("your_metadata.csv")

# Add metadata
adata.obs = meta

# Save in H5AD format
adata.write("data/exp1_counts_for_scenic_CELLLINE.h5ad")
```

## Note

The data files are not included in this repository. Please prepare your data according to the format above before running the pipeline.
