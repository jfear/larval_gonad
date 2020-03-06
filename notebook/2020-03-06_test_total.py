# %%
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon


# %%
df = pd.read_feather("../output/seurat3-cluster-wf/tpm_by_cluster.feather")

# %%
df['log'] = np.log1p(df.TPM)

# %%
res = []
for clus1, clus2  in combinations(df.cluster.cat.categories, 2):
    dat1 = df.query(f"cluster == '{clus1}'")['log']
    dat2 = df.query(f"cluster == '{clus2}'")['log']
    _, pval = mannwhitneyu(dat1, dat2, alternative="greater")
    res.append((clus1, clus2, pval))
res

# %%
