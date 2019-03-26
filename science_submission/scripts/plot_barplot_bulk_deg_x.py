import pickle

import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.stats import run_chisq

fname = snakemake.input.deg
fbgn2chrom = snakemake.input.fbgn2chrom
oname = snakemake.output[0]


def main():
    chroms = pd.read_parquet(fbgn2chrom)
    df = get_plot_data(chroms)
    pvals = get_pvals(chroms)

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(figsize=(1, 2))

    df.plot(kind='bar', stacked=True, legend=False, width=.9, color=['C0', 'gray', 'r'], ax=ax)
    ax.set_xticklabels([
        l.get_text().replace('chr', '')
        for l in ax.get_xticklabels()
    ], rotation=0, fontsize=8)

    if pvals.loc[('testis', 'flag_sig'), 'chrX'] == 1:
        ax.text(0, df.loc['chrX', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')

    if pvals.loc[('testis', 'flag_sig'), 'chr4'] == 1:
        ax.text(5, df.loc['chr4', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')

    ax.margins(0)
    sns.despine(ax=ax, left=True)
    ax.set_ylabel('% Genes')

    fig.savefig(oname, bbox_inches='tight')


def get_plot_data(chroms):
    bulk_sig = get_bulk_sig()

    df = bulk_sig.join(chroms).groupby('chrom').bias.value_counts().unstack()
    df = df.div(df.sum(axis=1), axis='rows') * 100
    return df.loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4'], ['testis', 'None', 'ovary']]


def get_pvals(chroms):
    bulk_sig = get_bulk_sig()
    df = (
        bulk_sig
        .join(chroms)
        .groupby('chrom').bias.value_counts()
        .unstack()
        .loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']]
        .T
    )
    return run_chisq(df)


def get_bulk_sig():
    bulk_sig = (
        pd.read_csv(fname, sep='\t', index_col=0)
            .dropna()
            .assign(testis_bias=lambda df: (df.log2FoldChange >= 1) & (df.padj <= 0.01))
            .assign(ovary_bias=lambda df: (df.log2FoldChange <= -1) & (df.padj <= 0.01))
    )

    bulk_sig.loc[bulk_sig.testis_bias, 'bias'] = 'testis'
    bulk_sig.loc[bulk_sig.ovary_bias, 'bias'] = 'ovary'
    bulk_sig.bias = bulk_sig.bias.fillna('None')
    return bulk_sig


if __name__ == '__main__':
    main()
