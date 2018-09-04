"""Figure 1. Identification, annotation and validation of clusters."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import seaborn as sns

from larval_gonad.plotting import flip_ticks
from plot_barchart_x2a_ratio import plot_barchart_x2a_ratio
from plot_boxplot_all_genes_by_chrom import plot_boxplot_all_genes_by_chrom
from plot_boxplot_common_genes_by_chrom import plot_boxplot_common_genes_by_chrom
from plot_boxplot_testis_bias_genes_by_chrom import plot_boxplot_testis_bias_genes_by_chrom
from plot_barchart_deg import plot_barchart_deg

boxplot_kws = dict(linewidth=.5)
boxplot_stripplot_kws = dict(linewidth=.5, size=3)
boxplot_line_kws = dict(linewidth=.5)
boxplot_defaults = dict(box_kws=boxplot_kws, line_kws=boxplot_line_kws, stripplot_kws=boxplot_stripplot_kws)

barplot_defaults = dict(linewidth=.5)

def tweak_boxplots(axes, x=False, title=False, row_label=None):
    for ax in axes:
        if x:
            ax.set_xticklabels([])
            ax.set_xlabel('')
            ax.xaxis.set_visible(False)

        if title:
            ax.set_title('')

        if row_label:
            ax = axes[-1]
            ax.text(1.1, 0.5, row_label, rotation=-90, ha='left', va='center', fontsize=10, transform=ax.transAxes)


def main(fig):
    # Make large grid
    gs = GridSpec(2, 1, height_ratios=[.1, 1], hspace=.25)
    gsCol = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[1, 0], height_ratios=[1, 1, 1, .5], hspace=0.03)

    # Left column plot axes
    axXA = fig.add_subplot(gs[0, 0])
    axXA.set_ylabel('X/A')

    # Add plots
    plot_barchart_x2a_ratio(axXA)
    All = plot_boxplot_all_genes_by_chrom(gsCol[0, 0], **boxplot_defaults)
    Common = plot_boxplot_common_genes_by_chrom(gsCol[1, 0], **boxplot_defaults)
    Bias = plot_boxplot_testis_bias_genes_by_chrom(gsCol[2, 0], **boxplot_defaults)
    Deg = plot_barchart_deg(gsCol[3, 0], bar_kws=barplot_defaults)

    # Tweak plots
    tweak_boxplots(All, x=True, row_label='All')
    tweak_boxplots(Common, x=True, title=True, row_label='Common')
    tweak_boxplots(Bias, x=True, title=True, row_label='Testis Bias')
    tweak_boxplots(Deg, title=True, row_label='DEG')

    # Add labels
    txt_defaults = dict(transform=fig.transFigure, fontweight='bold')
    plt.text(0.02, 0.9, 'A', **txt_defaults)
    plt.text(0.02, 0.75, 'B', **txt_defaults)
    plt.text(0.02, 0.56, 'C', **txt_defaults)
    plt.text(0.02, 0.37, 'D', **txt_defaults)
    plt.text(0.02, 0.20, 'E', **txt_defaults)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    main(fig)
    #fig.savefig(snakemake.output[0], bbox_inches='tight')
    fig.savefig(snakemake.output[0])
