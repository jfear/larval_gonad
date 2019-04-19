import pickle
import matplotlib

matplotlib.use('Agg')

import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.plotting import format_pval

gonia = snakemake.input.gonia
early = snakemake.input.early
mid = snakemake.input.mid
late = snakemake.input.late
fbgn2chrom = snakemake.input.fbgn2chrom
background = snakemake.input.background

oname = snakemake.output[0]

GERMLINE = ['SP', 'E1°', 'M1°']
QUERY_MAPPER = {
    'SP': {
        'fname': gonia,
        'query': 'p_val_adj <= 0.01 & avg_logFC > 0',
        'title': 'SP-Biased',
    },
    'E1°': {
        'fname': early,
        'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        'title': 'E1°-Biased',
    },
    'M1°': {
        'fname': mid,
        'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        'title': 'M1°-Biased',
    },
}


def main():
    chroms = get_chroms()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, axes = plt.subplots(1, 3, figsize=(4, 1), sharex=True, sharey=True, gridspec_kw=dict(hspace=.18))

    for name, ax in zip(GERMLINE, axes):
        # Get DEG
        significant_fbgns = get_significant_fbgns(name)
        ct = make_crosstab(significant_fbgns, chroms)

        # Plot proportion of genes
        proportion_genes_significant = ct.div(ct.sum()).loc[True, ['X', 'A', '4']]
        proportion_genes_significant.plot(kind='bar', width=.8, edgecolor='k', lw=.5,
                                          color=['darkgray', 'w', 'darkgray'], ax=ax)

        # Add * for significance
        _, xpval = fisher_exact(ct[['X', 'A']])
        _, fourpval = fisher_exact(ct[['4', 'A']])
        format_pval(ax, 0, proportion_genes_significant['X'], xpval, fontsize=7)
        format_pval(ax, 2, proportion_genes_significant['4'], fourpval, fontsize=7)

        # Clean up axes
        sns.despine(ax=ax)
        ax.set(ylim=(0, .25), ylabel='Proportion of Genes', xlabel='')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

        # Add title
        ax.text(0.5, 1, QUERY_MAPPER[name]['title'], transform=ax.transAxes, ha='center', va='top', fontsize=6)

    fig.savefig(oname, bbox_inches='tight')


def get_chroms():
    with open(background, 'rb') as fh:
        bg = pickle.load(fh)

    return (
        pd.read_parquet(fbgn2chrom)
        .replace({
            'chrX': 'X',
            'chr2L': 'A',
            'chr2R': 'A',
            'chr3L': 'A',
            'chr3R': 'A',
            'chr4': '4',
        })
        .reindex(bg)
        .query('chrom != ["chrM", "chrY"]')
    )


def get_significant_fbgns(name):
    return (
        pd.read_csv(QUERY_MAPPER[name]['fname'], sep='\t', index_col=0)
        .rename_axis('FBgn')
        .query(QUERY_MAPPER[name]['query'])
        .index.tolist()
    )


def make_crosstab(fbgns, chroms):
    _chroms = chroms.copy()
    _chroms['biased'] = False
    _chroms.loc[_chroms.index.isin(fbgns), 'biased'] = True

    return (
        _chroms
        .groupby('biased').chrom.value_counts()
        .unstack()
        .fillna(0)
    )


if __name__ == '__main__':
    main()
