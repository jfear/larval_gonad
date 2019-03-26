import pickle

import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

fname = snakemake.input.deg
oname = snakemake.output[0]


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(1, 1, figsize=(4, 2))

    defaults = dict(alpha=.5, s=6, rasterized=True)

    # Plot non-sig
    df[~df.testis_bias & ~df.ovary_bias].plot('baseMean', 'log2FoldChange', kind='scatter', color='gray', ax=ax,
                                              zorder=0, **defaults)
    ax.text(.1, 10,
            f'No-Bias = {(~df.ovary_bias & ~df.testis_bias).sum():,} ({( ~df.ovary_bias & ~df.testis_bias).mean() * 100:.0f}%)',
            color='gray', fontsize=8)

    # Plot testis bias
    df[df.testis_bias].plot('baseMean', 'log2FoldChange', kind='scatter', ax=ax, zorder=10, **defaults)
    ax.text(.5, 15, f'Testis-Bias = {df.testis_bias.sum():,} ({df.testis_bias.mean() * 100:.0f}%)', color='C0',
            fontsize=8)

    # Plot ovary bias
    df[df.ovary_bias].plot('baseMean', 'log2FoldChange', kind='scatter', ax=ax, color='r', zorder=10, **defaults)
    ax.text(.5, -12, f'Ovary-Bias = {df.ovary_bias.sum():,} ({df.ovary_bias.mean() * 100:.0f}%)', color='r',
            fontsize=8)

    ax.set_xscale('log')
    sns.despine(ax=ax)
    ax.axhline(0, color='k', lw=1, ls='--')

    fig.savefig(oname, bbox_inches='tight')


def get_data():
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
