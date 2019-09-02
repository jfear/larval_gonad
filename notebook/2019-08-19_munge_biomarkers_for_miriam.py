#%%
import numpy as np
import pandas as pd

from more_itertools import flatten
from larval_gonad.config import read_config

import os
try:
    os.chdir(os.path.join(os.getcwd(), 'notebook'))
    print(os.getcwd())
except:
    pass


#%%
# import gene annotations and lit genes
config = read_config('../config/common.yaml')
lit_genes = list(flatten(read_config('../config/literature_genes.yaml').values()))
fbgn2symbol = pd.read_feather("../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()

#%%
# import biomarkers and filter for alpha <= 0.01
biomarkers = (
    pd.read_feather("../output/seurat3-cluster-wf/combined_n3_biomarkers.feather")
    .assign(cluster=lambda x: x.cluster.map(config['cluster_annot']))
    .query("p_val_adj <= 0.01")
)

#%%
# Unstack to flag if a gene was a biomaker for a cluster
wide = pd.pivot_table(biomarkers, values="avg_logFC", columns="cluster", index='FBgn')
wide[wide.notnull()] = True
wide[wide.isnull()] = False

#%%
# Re-order columns and keep only lit genes
df_out = (
    wide.loc[:, config['cluster_order']]
    .reindex(lit_genes).fillna(0)
    .astype(int)
    .join(fbgn2symbol)
    .set_index("gene_symbol", append=True)
)

#%%
df_out.to_excel("../output/notebook/2019-08-19_munge_biomarkers_for_miriam.xlsx")




#%%
