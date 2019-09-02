import matplotlib

matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

FNAME = snakemake.input.tpm
METADATA = snakemake.input.metadata
GENE_METADATA = snakemake.input.gene_metadata

CLUSTER_COLORS = snakemake.params.cluster_colors

ONAME = snakemake.output[0]

# Debug Settings
# FNAME = 'output/science_submission/tpm_by_cluster_rep.feather'
# METADATA = 'output/seurat3-cluster-wf/combined_n3_metadata.feather'
# GENE_METADATA = 'references/gene_annotation_dmel_r6-26.feather'
# import yaml
# config = yaml.safe_load(open('config/common.yaml'))
# CLUSTER_COLORS = yaml.full_load(open('config/colors.yaml'))['clusters']


def main():
    y_fbgns = pd.read_feather(GENE_METADATA).query('FB_chrom == "Y"').FBgn.values
    tpm = pd.read_feather(FNAME).pipe(lambda df: df[df.FBgn.isin(y_fbgns)])

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots()
    sns.boxplot('cluster', 'TPM', data=tpm, 
        palette=CLUSTER_COLORS, 
        ax=ax, 
        linewidth=.5, 
        showfliers=False,
        notch=True,
    )
    ax.set(xlabel='', ylabel='Normalized Count (TPM)')

    fig.savefig(ONAME, bbox_inches='tight')


if __name__ == '__main__':
    main()
