import scanpy as sc
import numpy as np
import pandas as pd
from scipy import io
import pertpy as pt
from pathlib import Path

adata = sc.read_h5ad(Path(PATH_FILE)/'OSMGT_v1.h5ad')
# add pca
from sklearn.decomposition import PCA
pca = PCA(n_components=20)
pca = pca.fit_transform(adata.obsm['rss'])
adata.obsm['rss_pca'] = pca

# add cell line
metadata = pd.read_csv(Path(PATH_FILE)/"meta_OSMGT.tsv", sep='\t')
adata.obs['Cell_Line'] = metadata.Cell_Line
adata.obs['cell_condition'] = adata.obs['Cell_Line'].astype(str)+'_'+  adata.obs['Condition'].astype(str)

# Define distances
distance = pt.tl.Distance("edistance", obsm_key="rss_pca")

distance_mmd = pt.tl.Distance("mmd", obsm_key="rss_pca")

# Calculate distances
df = distance_mmd.pairwise(adata, groupby='cell_condition',n_jobs=15)
df.to_csv(Path(PATH_FILE)/'/mmd_general_OSMGT_self_pca.tsv',sep='\t')
print('MMD done')


df = distance.pairwise(adata, groupby='cell_condition',n_jobs=15)
df.to_csv(Path(PATH_FILE)/'e_distance_general_OSMGT_self_pca.tsv',sep='\t')
print('E-dist done')