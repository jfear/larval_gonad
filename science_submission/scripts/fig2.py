"""Figure 1. Identification, annotation and validation of clusters."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from larval_gonad.plotting import flip_ticks

from plot_heatmap_X_genes import plot_heatmap_X_genes
from plot_heatmap_Y_genes import plot_heatmap_Y_genes
from plot_heatmap_X_genes_pseudotime import plot_heatmap_X_genes_pseudotime
from plot_heatmap_Y_genes_pseudotime import plot_heatmap_Y_genes_pseudotime
from plot_barchart_germ_x2a_ratio import plot_barchart_germ_x2a_ratio
from plot_scatter_x2a_pseudotime import plot_scatter_x2a_pseudotime


def main(fig):
    # Make large grid
    gs = GridSpec(2, 3, width_ratios=[1, 1, 1], wspace=0.2)

    # Make X and Y heatmaps section for Seurat Clusters
    gsH1 = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[0, 0], width_ratios=[0.1, 1, 1], height_ratios=[1, .1], hspace=0)

    axCbar1 = fig.add_subplot(gsH1[0, 0])

    axX = fig.add_subplot(gsH1[0, 1])
    axXLabel = fig.add_subplot(gsH1[1, 1])

    axY = fig.add_subplot(gsH1[0, 2])
    axYLabel = fig.add_subplot(gsH1[1, 2])

    # Make X and Y heatmaps section for Psuedotime
    gsH2 = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[1, 0], width_ratios=[0.1, 1, 1], height_ratios=[1, .1], hspace=0)

    axCbar2 = fig.add_subplot(gsH2[0, 0])

    axXPT = fig.add_subplot(gsH2[0, 1])
    axXPTLabel = fig.add_subplot(gsH2[1, 1])

    axYPT = fig.add_subplot(gsH2[0, 2])
    axYPTLabel = fig.add_subplot(gsH2[1, 2])

    # Make grid for X2A section
    gsX2A = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0, 1], height_ratios=[1, 1, .1])
    axX2A = fig.add_subplot(gsX2A[0, 0])

    gsPT = GridSpecFromSubplotSpec(2, 1, subplot_spec=gsX2A[1:, 0], height_ratios=[1, .2], hspace=0)
    axX2APT = fig.add_subplot(gsPT[0, 0])
    axX2APTLabel = fig.add_subplot(gsPT[1, 0])

    # Other sections
    axQ = fig.add_subplot(gs[1, 1])
    axDEG = fig.add_subplot(gs[0, 2])
    axGO = fig.add_subplot(gs[1, 2])

    # Add plots
    plot_heatmap_X_genes(axX, axXLabel, cbar_ax=axCbar1)
    plot_heatmap_Y_genes(axY, axYLabel, cbar=False, yticklabels=False)

    plot_heatmap_X_genes_pseudotime(axXPT, axXPTLabel, cbar_ax=axCbar2)
    plot_heatmap_Y_genes_pseudotime(axYPT, axYPTLabel, cbar=False, yticklabels=False)

    plot_barchart_germ_x2a_ratio(axX2A)
    plot_scatter_x2a_pseudotime(axX2APT, axX2APTLabel, s=2, alpha=.6)

    axQ.axis('off')
    axQ.text(0.5, 0.5, '??', ha='center', va='center')

    axDEG.axis('off')
    axDEG.text(0.5, 0.5, 'Differential Expression', ha='center', va='center')

    axGO.axis('off')
    axGO.text(0.5, 0.5, 'Functional Analysis', ha='center', va='center')

    # Tweak plots
    flip_ticks(axCbar1)
    flip_ticks(axCbar2)

#     flip_ticks(axY, pos='right')
#     plt.setp(axY.get_yticklabels(), rotation=0, fontsize=6)
#
#     flip_ticks(axYPT, pos='right')
#     plt.setp(axYPT.get_yticklabels(), rotation=0, fontsize=6)

    axX2A.set_xticklabels([])
    axX2A.set_xticks([])
    axX2APT.set_ylabel('')

    # Add labels
#     txt_defaults = dict(transform=fig.transFigure, fontweight='bold')
#     plt.text(0.09, 0.89, 'A', **txt_defaults)
#     plt.text(0.09, 0.50, 'B', **txt_defaults)
#     plt.text(0.43, 0.89, 'C', **txt_defaults)
#     plt.text(0.64, 0.89, 'D', **txt_defaults)
#     plt.text(0.64, 0.50, 'E', **txt_defaults)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8.5, 6))
    main(fig)
    fig.savefig(snakemake.output[0], bbox_inches='tight')
