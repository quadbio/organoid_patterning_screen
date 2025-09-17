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


query_parse = sc.read_h5ad('/links/groups/treutlein/USERS/nazbukina/fatima_revision/data/exp1_raw.h5ad')
query_parse.obs['medium_morphogen'] = query_parse.obs['medium'].astype(str) + '_' +  query_parse.obs['morphogen_full'].astype(str) 

out2 = "/links/groups/treutlein/USERS/nazbukina/fatima_revision/data/parse_2_HNOCA/"


mapper2 = hmap.AtlasMapper(hnoca_model)
mapper2.map_query(
    query_parse, 
    retrain='partial',
  max_epochs=100,
   
    batch_size=1024
)
mapper2.save(out2)

mapper2.compute_wknn(k=50)

X_latent = mapper2.get_latent_representation(query_parse)
query_parse.obsm["X_latent_HNOCA"] = X_latent

presence_scores = mapper2.get_presence_scores(split_by="medium_morphogen")


np.save(os.path.join(out2, "lat_rep_query.npy"), X_latent)
import pickle
with open(os.path.join(out2, "presence.pickle"), 'wb') as f:
    pickle.dump(obj=presence_scores, file=f)

mapper2.save(out2)