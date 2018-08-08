"""Figure 1. Identification, annotation and validation of clusters."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from larval_gonad.plotting import flip_ticks

from plot_tsne import plot_tsne
from plot_heatmap_all_genes import plot_heatmap_all_genes
from plot_heatmap_lit_genes import plot_heatmap_lit_genes
from plot_heatmap_ptrap_genes import plot_heatmap_ptrap_genes


def main(fig):
    # Make large grid
    gs = GridSpec(4, 3, width_ratios=[1, 1, 1], height_ratios=[1, .1, 1, .1], hspace=0, wspace=0.1)

    # Split colorbars off of the all heatmap grid so I can control spacing
    gs2 = GridSpecFromSubplotSpec(4, 3, subplot_spec=gs[:, 1], width_ratios=[.3, .1, 1],
                                  height_ratios=[1, .1, 1, .1], hspace=0, wspace=0.05)

    axtSNE = fig.add_subplot(gs[2, 0])

    axCbar1 = fig.add_subplot(gs2[0, 1])
    axCbar2 = fig.add_subplot(gs2[2, 1])

    axAll = fig.add_subplot(gs2[:3, 2])
    axAllLabel = fig.add_subplot(gs2[3, 2])

    axLit = fig.add_subplot(gs[0, 2])
    axLitLabel = fig.add_subplot(gs[1, 2])

    gsPtrap = gs[2, 2]
    axPtrapLabel = fig.add_subplot(gs[3, 2])

    # Add plots
    plot_tsne(axtSNE)
    plot_heatmap_all_genes(axAll, axAllLabel, cbar_ax=axCbar1)
    plot_heatmap_lit_genes(axLit, axLitLabel, cbar=False)
    plot_heatmap_ptrap_genes(gsPtrap, axPtrapLabel, expression_heatmap_kws=dict(cbar=False),
                             ptrap_heatmap_kws=dict(cbar_ax=axCbar2),
                             )

    # Tweak plots
    flip_ticks(axCbar1)
    flip_ticks(axCbar2)
    flip_ticks(axLit, pos='right')
    plt.setp(axLit.get_yticklabels(), rotation=0)

    # Add labels
    txt_defaults = dict(transform=fig.transFigure, fontweight='bold')
    plt.text(0.09, 0.89, 'A', **txt_defaults)
    plt.text(0.09, 0.50, 'B', **txt_defaults)
    plt.text(0.43, 0.89, 'C', **txt_defaults)
    plt.text(0.64, 0.89, 'D', **txt_defaults)
    plt.text(0.64, 0.50, 'E', **txt_defaults)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8.5, 6))
    main(fig)
    fig.savefig(snakemake.output[0], bbox_inches='tight')
