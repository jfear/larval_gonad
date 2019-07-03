import matplotlib
matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

GENE_METADATA = snakemake.input.gene_metadata
EVIDENCE = snakemake.input.lit_evidence
CMAP = snakemake.params.cmap

ONAME = snakemake.output[0]

# Debug Settings
# import os
# os.chdir("science_submission/scripts")
# EVIDENCE = '../../data/external/lit_genes.tsv'
# GENE_METADATA = '../../references/gene_annotation_dmel_r6-24.feather'
# CMAP = 'viridis'

def main():
    fbgn2symbol = pd.read_feather(GENE_METADATA, columns=['FBgn', 'gene_symbol']).set_index('FBgn').to_dict()['gene_symbol']
    evidence = pd.read_csv(EVIDENCE, sep='\t', index_col=0).rename(fbgn2symbol)

    colors = sns.color_palette('viridis', n_colors=100)

    plt.style.use("scripts/paper_1c.mplstyle")
    fig, ax = plt.subplots(figsize=(4, 8))
    sns.heatmap(
        evidence.iloc[:, :-1], 
        cmap=[colors[0], 'lightgray', colors[-1]],
        yticklabels=True,
        xticklabels=True,
        cbar=False,
        square=True,
        ax=ax
    )
    plt.setp(ax.get_yticklabels(), fontsize=8, fontfamily='Helvetica')

    # Add legend
    low = mpatches.Patch(color=colors[0], label = 'Not Expressed')
    high = mpatches.Patch(color=colors[-1], label = 'Expressed')
    none = mpatches.Patch(color='lightgray', label = 'Not Tested')
    ax.legend(loc='upper left', bbox_to_anchor=[1, 1], handles=[low, high, none])

    # Clean up X axis
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_tick_params(pad=0, length=2)

    # Clean up Y axis
    ax.set_ylabel("")
    ax.yaxis.set_tick_params(pad=0.1, length=2)
    plt.setp(ax.get_yticklabels(), fontstyle="italic", fontfamily="Helvetica", rotation=0, va="center")

    fig.savefig(ONAME, bbox_inches="tight")


if __name__ == "__main__":
    main()