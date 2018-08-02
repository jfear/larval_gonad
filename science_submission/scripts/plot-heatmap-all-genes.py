"""Plot heatmap of all genes.

"""

import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels
import common

zscores = pd.read_parquet(snakemake.input[0])
oname = snakemake.output[0]

# calculate linkages
link = linkage(zscores.values, 'average')
tree = dendrogram(link, no_plot=True)
leaves = tree['leaves']

# plot
fig = plt.figure(figsize=(2, 4))
gs = GridSpec(2, 2, width_ratios=[.1, 1], height_ratios=[1, .08], hspace=0)
axCbar = plt.subplot(gs[0, 0])
axMain = plt.subplot(gs[0, 1])
axLabel = plt.subplot(gs[1, 1])

sns.heatmap(zscores.iloc[leaves], cmap=config['colors']['heatmap'], yticklabels=False, xticklabels=False,
            vmin=-3, vmax=3, ax=axMain, cbar_ax=axCbar, rasterized=True)

add_color_labels(axLabel)

axMain.set_ylabel('')
axCbar.yaxis.set_ticks_position('left')
axCbar.yaxis.set_label_position('left')
axCbar.set_ylabel('Normalized Read Counts (Z-score)')

fig.savefig(oname)
