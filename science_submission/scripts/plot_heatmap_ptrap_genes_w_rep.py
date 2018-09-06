"""Plot heatmap of protein trap genes.

Plots a heatmap of tpm zscores of genes that we have protein trap scores for.
Protein traps were selected based on their presence in the biomarkers
identified by Seurat. Protein traps were imaged and scored by multiple people
and were summarized by Miriam into a score ranging from 0-10. Currently I am
using 6 ptraps based on what images I had from miriam.
"""

import pandas as pd
from scipy import ndimage
import matplotlib as mpl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels_w_rep, flip_ticks
from common import fbgn2symbol

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_heatmap_ptrap_genes_w_rep(gsMain, label_size=5, expression_heatmap_kws=None, ptrap_heatmap_kws=None):
    zscores = pd.read_parquet('../scrnaseq-wf/data/tpm_zscore_w_rep.parquet', columns=config['sel_cluster_order_w_rep'])
    ptrap_scores = pd.read_parquet('data/ptrap_scores.parquet')
    ptrap_scores.index = ptrap_scores.index.droplevel(0)

    # Expand Cyst category to have early mid and late
    ptrap_scores = ptrap_scores.rename({'Cyst Cells': 'Early Cyst Cells'}, axis=1)
    ptrap_scores['Mid Cyst Cells'] = ptrap_scores['Early Cyst Cells']
    ptrap_scores['Late Cyst Cells'] = ptrap_scores['Early Cyst Cells']
    ptrap_scores = ptrap_scores[config['sel_cluster_order']]

    # Expand the ptrap scores to have three replicates for each column
    ptrap_scores.columns = [x + '-rep1' for x in ptrap_scores.columns]
    for col in ptrap_scores.columns:
        for rep in [2, 3]:
            name = col.replace('rep1', 'rep{}'.format(rep))
            ptrap_scores[name] = ptrap_scores[col]

    ptrap_scores = ptrap_scores[config['sel_cluster_order_w_rep']]

    # Get list of genes with ptrap scores
    ptrap_genes = ptrap_scores.index

    # Pull out ptrap genes that we have scores for
    zscores.index = zscores.index.map(fbgn2symbol)
    zscores = zscores.reindex(ptrap_genes)

    # plot
    defaults = dict(yticklabels=False, xticklabels=False, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    if isinstance(expression_heatmap_kws, dict):
        defaults.update(expression_heatmap_kws)

    # make set set of defaults for protein scores plot
    defaults2 = dict(yticklabels=False, xticklabels=False, vmin=0, vmax=8, rasterized=True,
                     cmap='inferno', cbar_kws=dict(label='Arbitrary Score'))

    if isinstance(ptrap_heatmap_kws, dict):
        defaults2.update(ptrap_heatmap_kws)

    # plot
    fig = plt.gcf()
    gs = GridSpecFromSubplotSpec(1, 3, subplot_spec=gsMain, wspace=0.03)

    fnames = {
        'SRPK': '../data/external/miriam/65332_srpk_Image 2_c1+2+3.tif',
        'bol': '../data/external/miriam/5.18.18_64431_bol_Image 2-Image Export-03_c1+2+3.tif',
        'Piezo': '../data/external/miriam/60209_piezo_Image 1_retake_c1+2+3.tif',
    }

    gss = [gs[0, 0], gs[0, 1], gs[0, 2]]
    for gene, gs0 in zip(['SRPK', 'bol', 'Piezo'], gss):
        gs00 = GridSpecFromSubplotSpec(4, 1,
                                       subplot_spec=gs0,
                                       height_ratios=[.1, .1, .1, 1],
                                       hspace=0,
                                       )

        axbg = fig.add_subplot(gs00[:, 0], facecolor='k')
        axbg.set_xticks([])
        axbg.set_yticks([])

        ax1 = fig.add_subplot(gs00[0, 0])
        ax2 = fig.add_subplot(gs00[1, 0])
        ax3 = fig.add_subplot(gs00[2, 0])
        ax4 = fig.add_subplot(gs00[3, 0])

        # plot the heatmaps and labels
        add_color_labels_w_rep(ax1, s=label_size)
        sns.heatmap(zscores.loc[gene].to_frame().T, ax=ax2, **defaults)
        sns.heatmap(ptrap_scores.loc[gene].to_frame().T, ax=ax3, **defaults2)

        # Plot miriams images
        img = plt.imread(fnames[gene])
        if gene == 'Piezo':
            img = ndimage.rotate(img, 90)

        ax4.imshow(img, aspect='equal')
        ax4.axis('off')
        ax4.text(0.01, -0.2, gene, transform=ax3.transAxes, color='#01ff07', ha='left', va='top',
                 fontsize=9, fontweight='bold')


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(4, 1, height_ratios=[1, 0.02, 0.04, 0.02], hspace=0.1)

    gsMain = gs[0, 0]
    axCbar1 = fig.add_subplot(gs[1, 0])
    axCbar2 = fig.add_subplot(gs[3, 0])

    plot_heatmap_ptrap_genes_w_rep(
        gsMain,
        label_size=5,
        expression_heatmap_kws=dict(cbar_ax=axCbar1, cbar_kws=dict(orientation='horizontal')),
        ptrap_heatmap_kws=dict(cbar_ax=axCbar2, cbar_kws=dict(orientation='horizontal')),
    )

    fig.savefig(snakemake.output[0], bbox_inches='tight')
