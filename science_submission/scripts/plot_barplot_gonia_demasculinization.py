import pickle

import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.stats import run_chisq
from larval_gonad.plotting import format_pval

bulk_fname = snakemake.input.bulk_deg
sc_fname = snakemake.input.sc_deg
fbgn2chrom = snakemake.input.fbgn2chrom
oname = snakemake.output[0]


def main():
    chroms = pd.read_parquet(fbgn2chrom)
    df = get_plot_data(chroms)

    # I am using two axes to simulate an axis break
    plt.style.use('scripts/figure_styles.mplstyle')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(1, 2), gridspec_kw=dict(hspace=0.1), sharex=True)

    # Plot proportion
    df.plot(kind='bar', stacked=True, legend=True, width=.9, color=['w', 'darkgray', 'k'], edgecolor='k', lw=.3, ax=ax1)
    df.plot(kind='bar', stacked=True, legend=False, width=.9, color=['w', 'darkgray', 'k'], edgecolor='k', lw=.3, ax=ax2)
    ax1.legend(['Testis-Biased', 'Not Biased', 'Ovary-Biased'], loc='upper left', bbox_to_anchor=[1, 1], frameon=False)

    # Add * for significance
    pvals = get_pvals(chroms)
    kwargs = dict(va='top', fontsize=7, color='k')
    format_pval(ax2, 0, df.loc['chrX', 'testis'] - .01, pvals.loc[('testis', 'fdr q-value'), 'chrX'], **kwargs)
    format_pval(ax2, 5, df.loc['chr4', 'testis'] - .01, pvals.loc[('testis', 'fdr q-value'), 'chr4'], **kwargs)

    # Add axis breaks
    d = 0.015
    kwargs = dict(transform=ax1.transAxes, clip_on=False, color='k', lw=1)
    ax1.plot((d, -d), (d, - d), **kwargs)
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((d, -d), (1 + d, 1 - d), **kwargs)

    # Clean up axes
    ax1.set(ylim=(.91, 1))
    ax1.xaxis.set_visible(False)
    ax1.margins(0)
    sns.despine(ax=ax1, bottom=True)

    ax2.set(ylim=(0, .09))
    ax2.set_xticklabels([l.get_text().replace('chr', '') for l in ax2.get_xticklabels()], rotation=0, fontsize=7)
    ax2.margins(0)
    sns.despine(ax=ax2)

    # Add y-label to span the two axis
    fig.text(-0.05, 0.5, 'Proportion of Genes', rotation=90, va='center', ha='right', fontsize=plt.rcParams['axes.labelsize'])
    ax1.set_title("Gonia-bias")

    fig.savefig(oname, bbox_inches='tight')


def get_plot_data(chroms):
    gonia_bias = get_gonia_biased()
    bulk_sig = get_bulk_sig().reindex(gonia_bias)

    df = bulk_sig.join(chroms, how='outer').fillna("None").groupby('chrom').bias.value_counts().unstack()
    df = df.div(df.sum(axis=1), axis='rows')
    return df.loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4'], ['testis', 'None', 'ovary']]


def get_pvals(chroms):
    gonia_bias = get_gonia_biased()
    bulk_sig = get_bulk_sig().reindex(gonia_bias)
    df = (
        bulk_sig
        .join(chroms, how='outer')
        .fillna("None")
        .groupby('chrom').bias.value_counts()
        .unstack()
        .loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']]
        .T
    )
    return run_chisq(df)


def get_bulk_sig():
    bulk_sig = (
        pd.read_csv(bulk_fname, sep='\t', index_col=0)
        .dropna()
        .assign(testis_bias=lambda df: (df.log2FoldChange >= 1) & (df.padj <= 0.01))
        .assign(ovary_bias=lambda df: (df.log2FoldChange <= -1) & (df.padj <= 0.01))
    )

    bulk_sig.loc[bulk_sig.testis_bias, 'bias'] = 'testis'
    bulk_sig.loc[bulk_sig.ovary_bias, 'bias'] = 'ovary'
    bulk_sig.bias = bulk_sig.bias.fillna('None')
    return bulk_sig


def get_gonia_biased():
    return (
        pd.read_csv(sc_fname, sep='\t', index_col=0)
        .query('p_val_adj <= 0.01 & avg_logFC > 0')
        .index.tolist()
    )


if __name__ == '__main__':
    main()
