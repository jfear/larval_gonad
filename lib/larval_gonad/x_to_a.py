"""Functions for use with X To A Analysis."""
import os
from pathlib import Path
from typing import Union, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns

from pysam import AlignmentFile

from .config import memory, config

CHROMS = ['X', '2L', '2R', '3L', '3R', '4', 'Y']
CHROMS_CHR = ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrY']

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


def x_autosome_boxplot(x, y, pvalue_cutoff=0.001, x_chrom=None, autosomes=None,
                       **kwargs):
    """Boxplot to compare X to Autsome"""
    if autosomes is None:
        autosomes = AUTOSOMES_CHR

    if x_chrom is None:
        x_chrom = 'chrX'

    _dat = kwargs.pop('data')
    chrX = _dat[_dat[x] == x_chrom][y]
    chrA = _dat[_dat[x].isin(autosomes)][y]
    stat, pvalue = mannwhitneyu(chrX, chrA, alternative='two-sided')

    med_x = chrX.median()
    med_a = chrA.median()
    diff = round(med_a - med_x, 2)

    # Pepare for plotting
    chrX = chrX.to_frame()
    chrX['class'] = 'X'

    chrA = chrA.to_frame()
    chrA['class'] = 'A'

    # Plot
    ax = sns.boxplot('class', y, data=pd.concat([chrX, chrA]),
                     order=['X', 'A'], **kwargs)
    ax.axhline(med_a, c=config['colors']['c2'], ls='--',
               label='Median Autsome Expression')
    ax.axhline(med_x, c=config['colors']['c3'], ls='--',
               label='Median X Expression')
    ax.legend()
    ax.text(.5, -.2, f'A - X = {diff}', transform=ax.transAxes, ha='center')

    # Add p-value
    if pvalue <= pvalue_cutoff:
        pvalue = clean_pvalue(pvalue)
        x1, x2 = 0, 1
        iqr = sns.utils.iqr(_dat[y])
        y, h, col = iqr + iqr * 2, .6, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*.5, y+h+.01, f"{pvalue}", ha='center',
                 va='bottom', color=col)

    return ax


def estimate_dcc(chrom_col: str,
                 count_col: str,
                 df: pd.DataFrame) -> Tuple[float, float, float]:
    """Estimate Dosage Compensation.

    Parameters
    ----------
    chrom_col: The column name with chromosomes info.
    count_col: The column name with the count info.
    df: Dataframe with chrom_col and count_col.

    Returns
    -------
    median of X, median of major autosomes, proporition compensation.

    """
    med_major = df.loc[df[chrom_col].isin(MAJOR_ARMS_CHR), count_col].median()
    med_x = df.loc[df[chrom_col] == 'chrX', count_col].median()

    try:
        prop_dcc = np.round((med_x / med_major) * 2, 2)
    except ZeroDivisionError:
        prop_dcc = np.nan

    return med_x, med_major, prop_dcc


def multi_chrom_boxplot(x, y, **kwargs):
    """Designed to be used with Seaborn.FacetGrid"""
    _dat = kwargs['data']
    _dat = _dat[_dat[x].isin(CHROMS_CHR)]
    curr_cluster = _dat['cluster'].values[0]
    num_cells = kwargs.pop('num_cells')[curr_cluster]
    num_genes = _dat.shape[0]

    ax = sns.boxplot(x, y, order=CHROMS_CHR, **kwargs)
    med_x, med_major, prop_dcc = estimate_dcc(x, y, _dat)
    ax.axhline(med_major, ls='--', color=config['colors']['c2'])

    # Clean up the pvalue for plotting
    pvalues = {}
    iqr = 0
    chromX = _dat[_dat.chrom == 'chrX']
    for g, df in _dat.groupby('chrom'):
        _iqr = sns.utils.iqr(df[y])
        if _iqr > iqr:
            iqr = _iqr
        if g == 'chrX':
            continue
        _, pval = mannwhitneyu(chromX[y], df[y], alternative='two-sided')
        if pval <= 0.001:
            pvalues[g] = pval

    multiplier = 2
    xloc = CHROMS_CHR.index('chrX')
    for k, v in pvalues.items():
        oloc = CHROMS_CHR.index(k)
        pval = clean_pvalue(v)
        y, h, col = iqr + iqr * multiplier, .4, 'k'
        plt.plot([xloc, xloc, oloc, oloc], [y, y+h, y+h, y], lw=1, c=col)
        plt.text((xloc+oloc)*.5, y+h+.01, f"{pval}", ha='center',
                 va='bottom', color=col)
        multiplier += 1

    ax.text(
        .05, .85,
        (f'X Compensation: {prop_dcc}\n'
         f'(n={num_cells} cells)\n'
         f'(g={num_genes:,} genes)'
         ),
        transform=ax.transAxes
    )

    # Add gene counts per chrom
#     gdsc = _dat.groupby(x)[y].describe()
#     gcnts = gdsc['count'].to_dict()
#     loc = (gdsc['75%'] + _dat.groupby(x)[y].apply(sns.utils.iqr) * 1.5).max()

#     for i, chrom in enumerate(CHROMS_CHR):
#         try:
#             pop_genes = np.round(gcnts[chrom] / num_genes * 100, 2)
#             ax.text(
#                 i,
#                 loc + .1, f'{gcnts[chrom]:0.0f}\n{pop_genes}%',
#                 ha='center'
#             )
#         except KeyError:
#             pass


def commonly_expressed(df, read_cutoff=0):
    """Common expressed genes.

    Determine if a gene is expressed above the read cutoff in more than 1/3 of
    cells.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame where rows are genes and columns are samples.
    read_cutoff : int
        Require more than cutoff reads. [0 default]

    Returns
    -------
    list :
        List generated by df.index.tolist() of genes that meet the expression
        criteria.

    """
    mask = (df > read_cutoff).sum(axis=1) > (df.shape[1] / 3)
    return df.index[mask].tolist()


def _get_mapper():
    REF_DIR = os.environ['REFERENCES_DIR']
    assembly = config['assembly']
    tag = config['tag']
    fname = Path(REF_DIR, assembly, tag, 'fb_annotation',
                 f'dmel_{tag}.fb_annotation')

    annot = pd.read_csv(fname, sep='\t')
    mapper = {}
    for i, row in annot.iterrows():
        fbgn = row.primary_FBgn
        mapper[fbgn] = fbgn

        if isinstance(row.annotation_ID, str):
            mapper[row.annotation_ID] = fbgn

        if isinstance(row.secondary_FBgn, str):
            for sec in row.secondary_FBgn.split(','):
                mapper[sec] = fbgn

        if isinstance(row.secondary_annotation_ID, str):
            for acn in row.secondary_annotation_ID.split(','):
                mapper[acn] = fbgn

    return mapper


MAPPER = _get_mapper()


def cleanup_FBgn(dat):
    """Convert FBgn or Accn to current FBgn."""
    res = []
    for k in dat:
        try:
            res.append(MAPPER[k])
        except KeyError:
            if isinstance(k, str):
                print(f'{k} not found in current annotation.')
    return np.array(res)


def get_gene_sets():
    PROJ_DIR = '../'

    genes = {}
    # Tau
    tau = pd.read_csv(Path(PROJ_DIR,
                           'output/2018-02-05_tau_haiwang_male_tau.tsv'),
                      sep='\t', index_col='FBgn', squeeze=True).dropna()
    tau_genes = cleanup_FBgn(tau[tau <= 0.4].index)
    genes['Haiwang_male_tau'] = tau_genes

    # TSPS
    tsps = pd.read_csv(Path(PROJ_DIR,
                            'output/2018-02-05_tau_haiwang_male_tsps.tsv'),
                       sep='\t', index_col='FBgn', squeeze=True).dropna()
    tsps_genes = cleanup_FBgn(tsps[tsps < 1.0].index)
    genes['Haiwang_male_tsps'] = tsps_genes

    # Naieve Bayes
    with open(Path(PROJ_DIR,
                   'data/external/Ferrari_et_al_2006_housekeeping_FBgn.txt')) as fh:
        housekeeping = np.array(fh.read().splitlines())

    bayes_genes = cleanup_FBgn(housekeeping)
    genes['Ferrari_housekeeping'] = bayes_genes

    # Protein groups from DroID and DPiM
    protein = pd.read_csv(Path(PROJ_DIR,
                               'data/external/DroID_DPiM_2018-03-29.txt'),
                          sep='\t', low_memory=False)

    network = defaultdict(set)
    for i, (bait, inter) in protein[['FBGN_BAIT', 'FBGN_INTERACTOR']].iterrows():
        network[bait].add(inter)
        network[inter].add(bait)

    # Sanity check to make sure I only have genes with at least one interaction
    dpim = set()
    for k, v in network.items():
        if len(v) > 1:
            dpim.add(k)
            dpim.union(set(v))

    dpim_genes = cleanup_FBgn(dpim)
    genes['DPiM_protein_complex'] = dpim_genes

    # Parse FlyBase Gene Groups
    header = [
        'FB_group_id',
        'FB_group_symbol',
        'FB_group_name',
        'Parent_FB_group_id',
        'Parent_FB_group_symbol',
        'FBgn',
        'gene_symbol',
    ]

    # Genes in any gene group
    fb_groups = pd.read_csv(Path(PROJ_DIR,
                                 'data/external/gene_group_data_fb_2017_03.tsv'),
                            sep='\t', comment='#', names=header)
    genes['flybase_groups'] = {}
    for i, group in fb_groups.groupby('FB_group_name'):
        fb_genes = cleanup_FBgn(group.FBgn)
        genes['flybase_groups'][i] = fb_genes

    return genes
