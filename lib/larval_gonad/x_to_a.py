"""Functions for use with X To A Analysis."""
from collections import defaultdict

import pandas as pd
from pysam import AlignmentFile

from .io import memory, config

CHROMS = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
CHROMS_CHR = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrY']

AUTOSOMES = ['2L', '2R', '3L', '3R', '4']
AUTOSOMES_CHR = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']

MAJOR_ARMS = ['2L', '2R', '3L', '3R']
MAJOR_ARMS_CHR = ['chr2L', 'chr2R', 'chr3L', 'chr3R']


def clean_pvalue(pval):
    if pval < 0.0001:
        pvalue = 'P<0.0001'
    else:
        pvalue = f'P={round(pval, 5)}'
    return pvalue


@memory.cache
def idx_stats_by_cluster(bam: str, cluster: str) -> pd.DataFrame:
    """Function to count the number of reads by cluster by chromosome.

    Parameters
    ----------
    bam : str
        The path to a bam file.
    cluster : str
        The path to a cluster file.

    Returns
    -------
    pd.DataFrame:
        DataFrame where rows are chromosome, columns are clusters, and values
        are number of aligned reads.

    """
    cluster = pd.read_csv(cluster, sep='\t')
    lookup = cluster.ident.to_dict()

    dat = AlignmentFile(bam, 'rb')
    results = defaultdict(lambda: defaultdict(int))

    for read in dat.fetch():
        if read.is_unmapped:
            continue

        chrom = read.reference_name
        if chrom not in CHROMS:
            continue

        try:
            cell = read.get_tag('CB').replace('-1', '')
            clus = lookup[cell]
            results[clus][chrom] += 1
        except KeyError:
            pass

    df = pd.DataFrame(results)
    df.index.name = 'chrom'
    return df
