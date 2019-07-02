import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns

GVCO = "/home/fearjm/Downloads/neo-X_pseudo/gonia_vs_cytes(1).tsv"
GVC = "../../output/seurat3-cluster-wf/combined_n3_gonia_vs_cytes.feather"
MULLERE = "../../output/neox-wf/sturgill_2007_mullerE.feather"
MULLERD = "../../output/neox-wf/sturgill_2007_mullerD.feather"

import yaml
config = yaml.safe_load(open('../../config/common.yaml'))
CLUSTER_ANNOT = config['cluster_annot']
CLUSTER_ORDER = config['cluster_order']

gvc = (
    pd.read_csv(GVCO, sep="\t")
    .set_index("primary_FBgn")
    .rename_axis("FBgn")
    .loc[:, "pct.2"]
    .rename("pct_cytes")
)

df = (
    pd.read_feather('../../output/seurat3-cluster-wf/raw_by_cluster_rep.feather')
    .groupby(['FBgn', 'cluster'])
    .UMI.sum()
    .unstack()
    .loc[:, ['EPS', 'PS1', 'PS2', 'PS3']]
    .sum(axis=1)
    .rename('total_cytes')
)

cytes = (
    pd.read_feather('../../output/seurat3-cluster-wf/combined_n3_metadata.feather', columns=['cell_id', 'cluster'])
    .set_index('cell_id')
    .assign(cluster=lambda x: x.cluster.cat.rename_categories(CLUSTER_ANNOT))
    .assign(cluster=lambda x: x.cluster.cat.reorder_categories(CLUSTER_ORDER))
    .query('cluster == ["EPS", "PS1", "PS2", "PS3"]')
    .index.tolist()
)

(
    pd.read_feather('../../output/cellselection-wf/raw.feather', columns=cytes)
    .groupby(['FBgn', 'cluster'])
    .UMI.sum()
    .unstack()
    .loc[:, ['EPS', 'PS1', 'PS2', 'PS3']]
    .sum(axis=1)
    .rename('total_cytes')
)

# Muller D Movement Calls
md = (
    pd.read_feather(MULLERD)
    .rename(columns={"Gene fate": "muller", "D.mel": "FBgn"})
    .set_index("FBgn")
    .join(df, how="left")
    .dropna()
)

md.muller.value_counts()

stats(md, name='total_cytes')

# Muller E Movement Calls
me = (
    pd.read_feather(MULLERE)
    .rename(columns={"Fate": "muller", "D.mel": "FBgn"})
    .set_index("FBgn")
    .join(df, how="left")
    .dropna()
)

me.muller.value_counts()

stats(me, name='total_cytes')


# Plots to match
plot_boxplot(md, me, 'total_cytes')


def stats(df, name='pct_cytes'):
    return mannwhitneyu(
        df.query('muller == "Conserved" | muller == "Move On"').loc[:, name],
        df.query('muller == "Gene Death" | muller == "Move Off"').loc[:, name],
        alternative="two-sided",
    )


def plot_boxplot(md, me, name='pct_cytes'):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 3), sharey=True)
    sns.boxplot("muller", name, data=md, order=["Conserved", "Move On", "Move Off", "Gene Death"], ax=ax1, showfliers=False)
    sns.boxplot("muller", name, data=me, order=["Conserved", "Move On", "Gene Death"], ax=ax2, showfliers=False)
    ax1.set(title="Muller D", xlabel="", ylabel="")
    plt.setp(ax1.get_xticklabels(), rotation=45)
    ax2.set(title="Muller E", xlabel="", ylabel="")
    plt.setp(ax2.get_xticklabels(), rotation=45)





if __name__ == "__main__":
    main()
