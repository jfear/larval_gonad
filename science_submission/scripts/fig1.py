"""Figure 1. Identification, annotation and validation of clusters."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from larval_gonad.plotting import flip_ticks, cluster_cmap

from plot_tsne import plot_tsne
from plot_heatmap_kmeans_all_genes import plot_heatmap_kmeans_all_genes
from plot_heatmap_lit_genes import plot_heatmap_lit_genes
from plot_heatmap_ptrap_genes import plot_heatmap_ptrap_genes
from plot_barchart_bulk_corr import plot_barchart_bulk_corr


def main(fig):
    # Make large grid
    gs = GridSpec(2, 3, width_ratios=[.7, .7, 1], wspace=0.05, hspace=0.05, left=0.01)

    axDia = fig.add_subplot(gs[0, 0])
    axtSNE = fig.add_subplot(gs[1, 0])

    # Heatmap All Genes:
    gsAll = GridSpecFromSubplotSpec(2, 3,
                                    subplot_spec=gs[:, 1],
                                    width_ratios=[0.15, .1, 1],
                                    height_ratios=[.05, 1],
                                    hspace=0,
                                    wspace=0.05
                                    )

    axCbar1 = fig.add_subplot(gsAll[1:, 1])
    axAllLabel = fig.add_subplot(gsAll[0, 2])
    axAll = fig.add_subplot(gsAll[1, 2])


    # Right column
    gsRight = GridSpecFromSubplotSpec(
        2, 3,
        subplot_spec=gs[:, 2],
        height_ratios=[.5, 1],
        width_ratios=[.1, 1, 1],
        hspace=0.1,
        wspace=0.1,
    )

    # Heatmap Lit Genes:
    gsLit = GridSpecFromSubplotSpec(2, 1,
                                    subplot_spec=gsRight[0, 1],
                                    height_ratios=[.2, 1],
                                    hspace=0
                                    )

    axLitLabel = fig.add_subplot(gsLit[0, 0])
    axLit = fig.add_subplot(gsLit[1, 0])

    # Bulk RNA-Seq:
    gsBulk = GridSpecFromSubplotSpec(2, 1,
                                     subplot_spec=gsRight[0, 2],
                                     height_ratios=[.2, 1],
                                     hspace=0
                                     )

    # Add a ylabel to the bulk axes area
    axBulk = fig.add_subplot(gsBulk[1, 0])
    flip_ticks(axBulk, pos='right')
    axBulk.xaxis.set_visible(False)
    axBulk.set_yticks([])
    axBulk.set_yticklabels([])
    axBulk.set_ylabel('Bulk RNA-Seq Correlation', labelpad=25)

    # Heatmap Ptrap Genes:
    gsPtrap = GridSpecFromSubplotSpec(1, 2,
                                      subplot_spec=gsRight[1, :],
                                      width_ratios=[1, .05],
                                      hspace=0,
                                      wspace=0.05,
                                      )

    axCbar2 = fig.add_subplot(gsPtrap[0, 1])
    gsPtrap2 = gsPtrap[0, 0]

    # Add plots
    img = plt.imread('../data/external/larval_testis_diagram.png')
    axDia.imshow(img)
    axDia.axis('off')

    plot_tsne(axtSNE)
    plot_heatmap_kmeans_all_genes(axAll, axAllLabel, cbar_ax=axCbar1)
    plot_heatmap_lit_genes(axLit, axLitLabel, cbar=False)
    plot_heatmap_ptrap_genes(gsPtrap2, label_size=3, expression_heatmap_kws=dict(cbar=False),
                             ptrap_heatmap_kws=dict(cbar_ax=axCbar2),
                             )
    axBulk1, axBulk2 = plot_barchart_bulk_corr(gsBulk[1, 0])

    # Tweak plots
    flip_ticks(axCbar1)
    plt.setp(axLit.get_yticklabels(), rotation=0)


    axBulk1.set_xticklabels([])
    axBulk1.set_xticks([])
    flip_ticks(axBulk1, pos='right')
    flip_ticks(axBulk2, pos='right')
    axBulk2.set_xlabel('Cell Cluster')

    # labels to distinguish cell types
    fontdict = dict(weight='bold', size=9, color=cluster_cmap['Early Cyst Cells'])
    axAllLabel.text(0.55, 1, 'Cyst', va='bottom', ha='left', transform=axAllLabel.transAxes, fontdict=fontdict)
    axLitLabel.text(0.55, 1, 'Cyst', va='bottom', ha='left', transform=axLitLabel.transAxes, fontdict=fontdict)

    fontdict.update(dict(color=cluster_cmap['Spermatogonia']))
    axAllLabel.text(0, 1, 'Germline', va='bottom', ha='left', transform=axAllLabel.transAxes, fontdict=fontdict)
    axLitLabel.text(0, 1, 'Germline', va='bottom', ha='left', transform=axLitLabel.transAxes, fontdict=fontdict)

    # Add labels
    txt_defaults = dict(transform=fig.transFigure, fontweight='bold')
    plt.text(-0.02, 0.89, 'A', **txt_defaults)
    plt.text(-0.02, 0.50, 'B', **txt_defaults)
    plt.text(0.30, 0.89, 'C', **txt_defaults)
    plt.text(0.53, 0.89, 'D', **txt_defaults)
    plt.text(0.73, 0.89, 'E', **txt_defaults)
    plt.text(0.53, 0.60, 'F', **txt_defaults)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8.5, 6))
    main(fig)
    fig.savefig(snakemake.output[0], bbox_inches='tight')
