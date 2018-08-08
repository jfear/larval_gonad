"""Plot heatmap of literature genes."""

import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels
import common

zscores = pd.read_parquet(snakemake.input[0])
lit_genes = snakemake.params.lit_genes
oname = snakemake.output[0]

# Pull out lit genes
zscores.index = zscores.index.map(common.fbgn2symbol)
zscores = zscores.reindex(lit_genes)

# plot
fig = plt.figure(figsize=(2, 2))
gs = GridSpec(2, 1, height_ratios=[1, .08], hspace=0)
axMain = plt.subplot(gs[0, 0])
axLabel = plt.subplot(gs[1, 0])

sns.heatmap(zscores, cmap=config['colors']['heatmap'], yticklabels=True, xticklabels=False,
            vmin=-3, vmax=3, ax=axMain, cbar=False)

add_color_labels(axLabel)
axMain.yaxis.set_ticks_position('right')
axMain.yaxis.set_label_position('right')
axMain.set_ylabel('')
plt.setp(axMain.get_yticklabels(), rotation=0, fontsize=10)

fig.subplots_adjust(left=0.1, right=0.80, bottom=0.01, top=0.9)
fig.savefig(oname)
