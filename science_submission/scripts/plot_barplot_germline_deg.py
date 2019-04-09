import pickle
import matplotlib

matplotlib.use('Agg')

import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

gonia = snakemake.input.gonia
early = snakemake.input.early
mid = snakemake.input.mid
late = snakemake.input.late
fbgn2chrom = snakemake.input.fbgn2chrom
background = snakemake.input.background

oname = snakemake.output[0]

AUTOSOMES = ['chr2L', 'chr2R', 'chr3L', 'chr3R']
COLORS = ['#e50000', '#ffffff', '#ffffff', '#ffffff', '#ffffff', '#e50000']
GERMLINE = ['SP', 'E1°', 'M1°', 'L1°']
QUERY_MAPPER = {
    'SP': {
        'fname': gonia,
        'query': 'p_val_adj <= 0.01 & avg_logFC > 0',
        'title': 'SP vs\nE1° M1° L1°',
        'yaxis': 'Proportion\nGonia-Baised'
    },
    'E1°': {
        'fname': early,
        'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        'title': 'SP vs E1°',
        'yaxis': 'Proportion\nE1°-Baised'
    },
    'M1°': {
        'fname': mid,
        'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        'title': 'E1° vs M1°',
        'yaxis': 'Proportion\nM1°-Baised'
    },
    'L1°': {
        'fname': late,
        'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        'title': 'M1° vs L1°',
        'yaxis': 'Proportion\nL1°-Baised'
    }
}


def main():
    chroms = get_chroms()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, axes = plt.subplots(4, 1, figsize=(1, 4), sharex=True, sharey=True, gridspec_kw=dict(hspace=.18))

    for name, ax in zip(GERMLINE, axes):
        significant_fbgns = get_significant_fbgns(name)
        ct = make_crosstab(significant_fbgns, chroms)
        _, xpval = fisher_exact(ct[['chrX', 'autosome']])
        _, fourpval = fisher_exact(ct[['chr4', 'autosome']])
        proportion_genes_significant = ct.div(ct.sum()).loc[True, ['chrX', ] + AUTOSOMES + ['chr4', ]]
        proportion_genes_significant.plot(kind='bar', width=.8, edgecolor='k', lw=.5, color=COLORS, ax=ax)

        # Add * for significance
        if xpval < 0.05:
            ax.text(0, proportion_genes_significant['chrX'], '*', ha='center', va='center')

        if fourpval < 0.05:
            ax.text(5, proportion_genes_significant['chr4'], '*', ha='center', va='center')

        # Add title
        ax.text(0.5, .8, QUERY_MAPPER[name]['title'], transform=ax.transAxes, ha='center', va='top', fontsize=6)
        # Clean up axes
        sns.despine(ax=ax)
        ax.set_ylim(0, .25)

        ax.set_xticklabels([
            l.get_text().replace('chr', '')
            for l in ax.get_xticklabels()
        ], rotation=0)

        ax.set_xlabel('')
        ax.set_ylabel(QUERY_MAPPER[name]['yaxis'])

    fig.savefig(oname, bbox_inches='tight')


def get_chroms():
    with open(background, 'rb') as fh:
        bg = pickle.load(fh)

    return pd.read_parquet(fbgn2chrom).reindex(bg).query('chrom != ["chrM", "chrY"]')


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
            .assign(autosome=lambda df: df[AUTOSOMES].sum(axis=1))
            .fillna(0)
    )


if __name__ == '__main__':
    main()
