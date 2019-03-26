import pickle

import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from larval_gonad.gene_ontology import run_flyslim

fname = snakemake.input.deg
background = snakemake.input.background
oname = snakemake.output[0]


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(1, 1, figsize=(2, 6))

    # Main plot
    scale = .5
    sc = ax.scatter(df['g_prop'], df['loc'], c=df['g_adj_pval'], s=df['g_cnt'] * scale, vmin=0, vmax=.1, cmap='coolwarm_r',
                    label='')
    ax.scatter(df['m_prop'], df['loc'], c=df['m_adj_pval'], s=df['m_cnt'] * scale, vmin=0, vmax=.1, cmap='coolwarm_r',
                    label='')

    ax.axvline(0, ls='--', lw=.5, color='k')

    ax.set_axisbelow(True)
    ax.grid(True, axis='y', alpha=.2)
    sns.despine(ax=ax, left=True, bottom=True)
    ax.margins(0.02)

    # Add color bar and legend
    plt.colorbar(sc, orientation='horizontal', ticks=[0, .05, .1], label='Adjusted p-value', fraction=.05, pad=.07)

    for sizes in [10, 100, 300]:
        plt.scatter([], [], c='k', alpha=0.3, s=sizes * scale, label=f'{sizes:,.0f}')

    legend = plt.legend(scatterpoints=1, frameon=False, labelspacing=.6, title='Number of Genes', loc='upper left',
                        bbox_to_anchor=[1, 1], bbox_transform=ax.transAxes)
    plt.setp(legend.get_title(), fontsize=8)

    # Clean up X
    ax.set_xlabel('Proportion of Genes')
    ax.set_xlim(-.2, .2)

    # Clean up Y axis
    ax.set_yticks(range(1, df.shape[0] + 1))
    ax.set_yticklabels(df.index.tolist())

    # Add titles
    ax.text(0.25, 1, 'SP-Biased', transform=ax.transAxes, ha='center', va='bottom');
    ax.text(0.75, 1, 'M1°-Biased', transform=ax.transAxes, ha='center', va='bottom');

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    with open(background, 'rb') as fh:
        bg = pickle.load(fh)
    deg = (
        pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_mid.tsv', sep='\t', index_col=0)
            .rename_axis('FBgn')
            .query('p_val_adj < 0.01')
    )

    gonia = run_go(deg.query('avg_logFC > 0').index.tolist(), bg, 'g_')
    mid = run_go(deg.query('avg_logFC < 0').index.tolist(), bg, 'm_')

    return (
        pd.concat([gonia, mid], axis=1, sort=True)
        .assign(min_pval=lambda df: df[['g_adj_pval', 'm_adj_pval']].min(axis=1))
        .query('min_pval <= 0.01')
        .assign(g_prop=lambda df: df['g_prop'] * -1)
        .assign(m_cnt=lambda df: df['m_cnt'].fillna(0))
        .assign(m_prop=lambda df: df['m_prop'].fillna(0))
        .assign(m_adj_pval=lambda df: df['m_adj_pval'].fillna(1))
        .assign(loc=lambda df: range(1, df.shape[0] + 1))
        .sort_values(by='g_prop')
    )


def run_go(fbgns, bg, prefix):
    cols = [prefix + 'cnt', prefix + 'prop', prefix + 'adj_pval']
    df = pd.DataFrame(columns=cols)
    for x in run_flyslim(fbgns, bg, cutoff=.9):
        if x.name == 'molecular_function':
            continue
        if x.NS == 'MF':
            df = df.append(pd.Series([x.study_count, x.ratio_in_study[0] / x.ratio_in_study[1], x.p_fdr_bh], index=cols,
                                     name=x.name))
    return df


if __name__ == '__main__':
    main()
