import os
import scarches
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import hnoca.map as hmap
import hnoca.stats as stats





hnoca_adata = sc.read("/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_atlas/data/phase3_final_0516/ZH_processed_cleanup_highvar_counts.h5ad")
import torch
model_path = "/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_atlas/data/VAE_models/scpoli_hierarchical123"
model = scarches.models.scPoli
hnoca_model = model.load(model_path, hnoca_adata, map_location=torch.device('cpu'))


query_adata = sc.read_h5ad('/links/groups/treutlein/USERS/nazbukina/fatima_revision/data/OSMGT_raw.h5ad')
#query_adata.obs["batch"] = query_adata.obs["sample"].astype(str).copy()
metadata = pd.read_csv("/links/groups/treutlein/USERS/nazbukina/fatima_revision/data/meta_OSMGT.tsv", sep='\t')
query_adata.obs['merged_cluster_annotation'] = metadata.merged_cluster_annotation
query_adata.obs['Condition'] = metadata.Condition


out = "/links/groups/treutlein/USERS/nazbukina/fatima_revision/data/OSMGT_2_HNOCA/"



mapper = hmap.AtlasMapper(hnoca_model)
mapper.map_query(
    query_adata, 
    retrain='partial',
  max_epochs=100,
   
    batch_size=1024
)
mapper.save(out)

mapper.compute_wknn(k=50)

X_latent = mapper.get_latent_representation(query_adata)
query_adata.obsm["X_latent"] = X_latent

presence_scores = mapper.get_presence_scores(split_by="Condition")


np.save(os.path.join(out, "lat_rep_query.npy"), X_latent)
import pickle
with open(os.path.join(out, "presence.pickle"), 'wb') as f:
    pickle.dump(obj=presence_scores, file=f)
mapper.save(out)