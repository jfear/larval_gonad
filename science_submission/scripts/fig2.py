"""Figure 1. Identification, annotation and validation of clusters."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from larval_gonad.plotting import flip_ticks
from plot_boxplot_all_genes_by_chrom import plot_boxplot_all_genes_by_chrom
from plot_boxplot_common_genes_by_chrom import plot_boxplot_common_genes_by_chrom
from plot_boxplot_germ_common_genes_by_chrom import plot_boxplot_germ_common_genes_by_chrom
from plot_boxplot_upregulated_genes_by_chrom import plot_boxplot_upregulated_genes_by_chrom
from plot_manhatten import plot_manhatten
from plot_genome_view import plot_genome_view
from plot_genome_view_log import plot_genome_view_log


def tweak_boxplots(axes, x=False, title=False, row_label=None):
    for ax in axes:
        if x:
            ax.set_xticklabels([])
            ax.set_xlabel('')

        if title:
            ax.set_title('')

        if row_label:
            ax = axes[-1]
            ax.text(1.1, 0.5, row_label, rotation=-90, ha='left', va='center', fontsize=10, transform=ax.transAxes)


def main(fig):
    # Make large grid
    gs = GridSpec(2, 2, hspace=.2, wspace=.2)

    # Create layout
    gsXA = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[0, 0])
    gsView = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 1])
    ax1 = fig.add_subplot(gsView[0, 0])
    ax2 = fig.add_subplot(gsView[1, 0])

    # Add plots
    all_genes = plot_boxplot_all_genes_by_chrom(gsXA[0, 0])
    common_genes = plot_boxplot_common_genes_by_chrom(gsXA[1, 0])
    germ_common_genes = plot_boxplot_germ_common_genes_by_chrom(gsXA[2, 0])
    upregulated_genes = plot_boxplot_upregulated_genes_by_chrom(gsXA[3, 0])
    manhatten = plot_manhatten(gs[1, 0])

    plot_genome_view(ax1)
    plot_genome_view_log(ax2)

    # Tweak plots
    tweak_boxplots(all_genes, x=True, row_label='All')
    tweak_boxplots(common_genes, x=True, title=True, row_label='Common\nAll Cells')
    tweak_boxplots(germ_common_genes, x=True, title=True, row_label='Common\nGerm Cells')
    tweak_boxplots(upregulated_genes, title=True, row_label='New')

    # Add labels
    txt_defaults = dict(transform=fig.transFigure, fontweight='bold')
#     plt.text(0.09, 0.89, 'A', **txt_defaults)
#     plt.text(0.09, 0.50, 'B', **txt_defaults)
#     plt.text(0.43, 0.89, 'C', **txt_defaults)
#     plt.text(0.64, 0.89, 'D', **txt_defaults)
#     plt.text(0.64, 0.50, 'E', **txt_defaults)


if __name__ == '__main__':
    fig = plt.figure(figsize=(11, 8))
    main(fig)
    #fig.savefig(snakemake.output[0], bbox_inches='tight')
    fig.savefig(snakemake.output[0])
