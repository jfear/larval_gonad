"""Plot heatmap of unique biomarkers"""
import matplotlib

matplotlib.use("Agg")

from itertools import chain
from more_itertools import flatten
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, dendrogram

from larval_gonad.io import feather_to_cluster_rep_matrix

FNAME = snakemake.input.zscores
GENE_METADATA = snakemake.input.gene_metadata
BIOMARKERS = snakemake.input.biomarkers

CLUSTER_ANNOT = snakemake.params.cluster_annot
CLUSTER_ORDER = snakemake.params.cluster_order
CMAP = snakemake.params.cmap

ONAME = snakemake.output[0]

# Debug settings
FNAME = "output/science_submission/zscore_by_cluster_rep.feather"
GENE_METADATA = "references/gene_annotation_dmel_r6-24.feather"
BIOMARKERS = "output/seurat3-cluster-wf/combined_n3_biomarkers.feather"
import yaml
config = yaml.safe_load(open('config/common.yaml'))
CLUSTER_ANNOT = config['cluster_annot']
CLUSTER_ORDER = config['cluster_order']
CMAP = "viridis"


def main():

fbgn2symbol = (
    pd.read_feather(GENE_METADATA, columns=["FBgn", "gene_symbol"])
    .set_index("FBgn")
    .to_dict()["gene_symbol"]
)

biomarkers = (
    pd.read_feather(BIOMARKERS, columns=['FBgn', 'cluster'])
    .assign(cluster=lambda df: df.cluster.cat.rename_categories(CLUSTER_ANNOT))
    .assign(cluster=lambda df: df.cluster.cat.reorder_categories(CLUSTER_ORDER))
    .drop_duplicates(subset="FBgn", keep=False)
)

zscores = feather_to_cluster_rep_matrix(FNAME).reindex(biomarkers.FBgn)

# cluster genes bases on expression
link = linkage(zscores.values, 'average')
tree = dendrogram(link, no_plot=True)
leaves = tree['leaves']
zscores = zscores.iloc[leaves, :]

plt.style.use('scripts/figure_styles.mplstyle')
fig = plt.figure(figsize=(4, 8))
gs = GridSpec(2, 1, height_ratios=[1, .02])
ax = fig.add_subplot(gs[0, 0])
cax = fig.add_subplot(gs[1, 0])
sns.heatmap(
    zscores,
    xticklabels=True,
    yticklabels=False,
    vmin=-3,
    vmax=3,
    rasterized=False,
    cmap=CMAP,
    ax=ax,
    cbar_ax=cax,
    cbar_kws=dict(label='Z-Score (TPM)', ticks=[-3, 0, 3], orientation='horizontal')
)

# Clean up X axis
ax.set_xlabel('')
ax.xaxis.set_ticks_position('top')
ax.set_xticklabels(list(chain.from_iterable([('', x, '') for x in cluster_order])), ha='center', va='bottom')

# Add additional x annotations
yloc = 0 - (df.shape[0] * .05)
pad = yloc * .1
ax.text(6, yloc + pad, 'Germline', ha='center', fontsize=6, color=cluster_colors[0], va='bottom')
ax.text(17, yloc + pad, 'Somatic\nCyst', ha='center', fontsize=6, color=cluster_colors[4], va='bottom')
ax.text(24, yloc + pad, 'Somatic\nOther', ha='center', fontsize=6, color=cluster_colors[8], va='bottom')
ax.text(31, yloc + pad, 'Unknown', ha='center', fontsize=6, color=cluster_colors[-1], va='bottom')
lines = [
    plt.Line2D([0, 12], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
    plt.Line2D([12, 21], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
    plt.Line2D([21, 24], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
    plt.Line2D([24, 27], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
    plt.Line2D([27, 35], [yloc, yloc], color=cluster_colors[-1], lw=1.5, clip_on=False),
]

for l in lines:
    ax.add_line(l)

# Add lines separating cell types
for i in range(1, 12):
    ax.axvline(i * 3, color='w', ls='--', lw=.5)

# Clean up Y axis
ax.set_ylabel('')

# Add lines separating biomarker groups
for loc in np.cumsum(unique_genes_per_cluster)[:-1]:
    ax.axhline(loc, color='w', ls='--', lw=.5)

# Add additional y annotations
loc = 0
for i, cnt in enumerate(unique_genes_per_cluster):
    if i == 10:
        x = 38
    else:
        x = 35
    ax.text(-.5, loc + (cnt / 2), cluster_order[i], ha='right', va='center', fontweight='bold', fontsize=5.5)
    ax.text(x, loc + (cnt / 2), f'n={cnt}', ha='left', va='center', fontweight='bold', fontsize=5.5)
    loc += cnt

fig.savefig(oname, bbox_inches='tight')



if __name__ == '__main__':
    main()
