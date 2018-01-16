#!/usr/bin/env python
"""This script makes tSNE plots for all genes."""
import os
import sys
from pathlib import Path
import string

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
sys.path.insert(0, '../lib')
from larval_gonad.plotting import TSNEPlot, add_styles
from larval_gonad.logging import logger

REF = os.environ['REFERENCES_DIR']
add_styles('../config/stylelib')

def sanitize_fname(fname):
    valid_chars = "-_.%s%s" % (string.ascii_letters, string.digits)
    return ''.join([x for x in fname if x in valid_chars])

def plot_gene(data, fbgn, symbol, output, **kwargs):
    symbol = sanitize_fname(symbol)
    fname = str(Path(output, f'{fbgn}_{symbol}.png'))
    if Path(fname).exists():
        return

    df = data[['tSNE_1', 'tSNE_2', fbgn]]

    with plt.style.context(['paper-wide', 'default']):
        fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1.3, 1]})
        TSNEPlot('tSNE_2', 'tSNE_1', data=df, hue=fbgn, s=10,
                 ax=ax1, title='Normalized Expression\n(Continuous)', **kwargs)

        TSNEPlot('tSNE_2', 'tSNE_1', data=df, hue=df[fbgn] > 0,
                 cmap={
                     '0': 'w',
                     '1': 'k',
                 }, s=10, ax=ax2, alpha=.6, edgecolor='k', title='Normalized Expression\n(Binary)', **kwargs)

        fig.suptitle(f'{symbol} ({fbgn})');
        plt.tight_layout(rect=[0, 0, .9, .9])
        plt.savefig(fname)
        plt.close()

if __name__ == '__main__':
    # gene annotations
    fbgn2symbol = pd.read_csv(str(Path(REF, 'dmel/r6-16/fb_annotation/dmel_r6-16.fb_annotation')), sep='\t',
                              usecols=['gene_symbol', 'primary_FBgn'], index_col='primary_FBgn').to_dict()['gene_symbol']

    # Colors
    colors = sns.color_palette('Reds')
    color2 = sns.color_palette('Greys')
    colors[0] = color2[0]

    # testes
    logger.info('Plotting Testes')
    FIGS = '../output/figures/testis_tsne'
    Path(FIGS).mkdir(exist_ok=True)

    DAT = '../output/testes_scRNAseq_pilot'
    tsne = pd.read_csv(Path(DAT, 'tsne.tsv'), sep='\t')
    norm = pd.read_csv(Path(DAT, 'normalized_read_counts.tsv'), sep='\t')
    data = tsne.join(norm.T)

    for fbgn in data.columns[2:]:
        symbol = fbgn2symbol[fbgn]
        plot_gene(data, fbgn, symbol, FIGS, palette=colors)

    # combined
    logger.info('Plotting Combined')
    FIGS = '../output/figures/combined_tsne'
    Path(FIGS).mkdir(exist_ok=True)

    DAT = '../output/combined_ovary_and_testis_scRNAseq_pilot'
    tsne = pd.read_csv(Path(DAT, 'tsne.tsv'), sep='\t')
    norm = pd.read_csv(Path(DAT, 'normalized_read_counts.tsv'), sep='\t')
    data = tsne.join(norm.T)

    for fbgn in data.columns[2:]:
        symbol = fbgn2symbol[fbgn]
        plot_gene(data, fbgn, symbol, FIGS, palette=colors)
