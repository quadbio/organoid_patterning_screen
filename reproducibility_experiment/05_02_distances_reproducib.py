import scanpy as sc
import numpy as np
import pandas as pd
from scipy import io
import pertpy as pt
from pathlib import Path


adata = sc.read_h5ad(Path(PATH_FILE)/'exp1_processed_5_1.h5ad')

from sklearn.decomposition import PCA
pca = PCA(n_components=20)
parse_pca = pca.fit_transform(adata.obsm['rss'])
adata.obsm['rss_pca'] = parse_pca


# Define distance
distance = pt.tl.Distance("edistance", obsm_key="rss_pca")
distance_mmd = pt.tl.Distance("mmd", obsm_key="rss_pca")

adata = adata[adata.obs.cell_line.isin(['H1', 'H9', 'WTC', 'WIBJ2'])].copy()

adata.obs['cell_condition_batch'] = adata.obs['cell_line'].astype(str)+'_'+  adata.obs['medium_morphogen'].astype(str)+'_'+  adata.obs['batch'].astype(str)

# Calculate distance
df = distance_mmd.pairwise(adata, groupby='cell_condition_batch',n_jobs=15)

df.to_csv(Path(PATH_FILE)/'mmd_general.tsv',sep='\t')
print('MMD done')


df = distance.pairwise(adata, groupby='cell_condition_batch',n_jobs=15)
df.to_csv(Path(PATH_FILE)/'e_distance_general.tsv',sep='\t')

print('E-dist done')

