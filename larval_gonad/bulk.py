"""Helper functions for working with bulk data.

We performed bulk RNA-seq and this is a set of helpers for dealing with this
data.

"""
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
from lcdblib.plotting import maPlot, PairGrid, corrfunc

from .cell_selection import filter_gene_counts_by_barcode

TESTIS_BULK = [
    'B5_TCP',
    'B6_TCP',
    'B7_TCP',
    'B8_TCP',
]


def read_bulk(path, filter=None, pattern='*/*.featurecounts.txt'):
    """Read in a folder of feature count data.

    Using the lcdb-wf, featurecounts are organized in a set of sub-folders for
    each sample. Given a path will read in the data and return a dataframe.
    Optionally a list of sample names can be given to filter by.

    Parameters
    ----------
    path : str
        Directory path to output from the lcdb-wf.
    filter : None | list
        List of sample names to include. Defaults to use the TCP libraries.
    pattern : str
        Glob pattern for finding the featurecounts files.

    Example
    -------
    >>> df = read_build('../bulk-rnaseq-wf/data/rnaseq_samples',
            filter=['B5_TCP', 'B6_TCP'])
    """
    bulk = Path(path)

    dfs = []
    for fname in bulk.glob(pattern):
        sname = fname.parent.name
        if (filter is not None) & (sname in filter):

            dat = pd.read_csv(fname, sep='\t', comment='#',
                              index_col=[0]).iloc[:, -1]

            dat.name = sname
            dfs.append(dat)

    bulk_dat = pd.concat(dfs, axis=1)
    bulk_dat = bulk_dat[bulk_dat.columns.sort_values()]
    return bulk_dat


def read_bulk_for_lengths(path, filter=None, pattern='*/*.featurecounts.txt'):
    """Read in a folder of feature count data to get gene lengths.

    Using the lcdb-wf, featurecounts are organized in a set of sub-folders for
    each sample. Given a path will read in the data and return a dataframe.
    Optionally a list of sample names can be given to filter by.

    Parameters
    ----------
    path : str
        Directory path to output from the lcdb-wf.
    filter : None | list
        List of sample names to include. Defaults to use the TCP libraries.
    pattern : str
        Glob pattern for finding the featurecounts files.

    Example
    -------
    >>> df = read_build('../bulk-rnaseq-wf/data/rnaseq_samples',
            filter=['B5_TCP', 'B6_TCP'])
    """
    bulk = Path(path)

    dfs = []
    for fname in bulk.glob(pattern):
        sname = fname.parent.name
        if (filter is not None) & (sname in filter):

            dat = pd.read_csv(fname, sep='\t', comment='#',
                              index_col=[0]).iloc[:, -2]

            dat.name = 'length'
            dfs.append(dat)

    bulk_dat = pd.concat(dfs, axis=0)
    return bulk_dat.to_frame().reset_index().drop_duplicates().set_index('Geneid').length


def plot_bulk_pairwise_corr(bulk_dat, subplots_kws=None, scatter_kws=None,
                            corrfunc_kws=None):
    """Plot a pairgrid of RNA-seq data.

    The upper triangle is the scatter plot and spearman correlation. The lower
    triangle is a common MA-Plot. The diagonal is the density.

    bulk_dat : pd.DataFrame
        DataFrame with RNA-seq data (genes, samples)

    """
    if subplots_kws is None:
        subplots_kws = {}

    if scatter_kws is None:
        scatter_kws = {}

    if corrfunc_kws is None:
        corrfunc_kws = {}

    subplots_default = {
        'sharex': False,
        'sharey': False
    }
    subplots_default.update(subplots_kws)

    scatter_default = {
        's': 10
    }
    scatter_default.update(scatter_kws)

    corrfunc_default = {
    }
    corrfunc_default.update(corrfunc_kws)

    g = PairGrid(bulk_dat, subplots_kws=subplots_default)
    g.map_lower(maPlot, scatter_kws=scatter_default)
    g.map_upper(plt.scatter, **scatter_default)
    g.map_upper(corrfunc, **corrfunc_default)
    g.map_diag(sns.kdeplot)

    return g


def scRNAseq_corr_distribution(umi, raw, bulk_dat, start=200,
                               interval=100, stop=10000):
    """Calculate the correlation distribution between scRNASeq and Bulk.

    Iterate by intervals of cells and calculate the correlation of summed
    scRNASeq vs Bulk RNA-Seq.

    Parameters
    ----------
    umi : pd.DataFrame
        DataFrame of UMI counts by Cell (tidy)
    raw : CellRangerCounts
        A named tuple of CellRangerCounts.
    bulk_dat : pd.DataFrame
        DataFrame of bulk RNA-seq data (genes, samples)
    start : int
        Number of cells to start with [default 200]
    interval : int
        Number of cells to add each iteration [default 100]
    stop : int
        Number of cells to stop at [default 10,000]

    Returns
    -------
    pd.DataFrame
        Rows are the number of UMI sorted cells. Columns are Bulk RNASeq
        samples. Values are Spearman r coefficients.

    """
    _umi = umi.sort_values(by='umi_count', ascending=False)

    res = []
    loc = start
    while loc < stop:
        dat = filter_gene_counts_by_barcode(_umi.index[:loc], raw).sum(axis=1)
        corrs = []
        for col in bulk_dat.columns:
            corrs.append(spearmanr(bulk_dat[col], dat).correlation)

        res.append([loc, *corrs])
        loc += interval

    col_names = ['Cell Number']
    col_names.extend(bulk_dat.columns)

    df = pd.DataFrame(res, columns=col_names)

    return df.set_index('Cell Number')


def plot_corr_distribution(corr):
    fig, axes = plt.subplots(2, 2, sharex=True)

    for col, ax in zip(corr.columns, axes.flatten()):
        ax.plot(corr[col])
        ax.set_title(col)
        ax.set_ylabel('Spearman r')
        ax.set_xlabel('Cells')

    plt.tight_layout()


def scRNAseq_corr_distribution_random(umi, raw, bulk_dat, interval=100,
                                      stop=10000, random_state=42):
    """Calculate the correlation distribution between scRNASeq and Bulk.

    Iterate by intervals of cells and calculate the correlation of summed
    scRNASeq vs Bulk RNA-Seq.

    Parameters
    ----------
    umi : pd.DataFrame
        DataFrame of UMI counts by Cell (tidy)
    raw : CellRangerCounts
        A named tuple of CellRangerCounts.
    bulk_dat : pd.DataFrame
        DataFrame of bulk RNA-seq data (genes, samples)
    interval : int
        Number of cells to add each iteration [default 100]
    stop : int
        Number of cells to stop at [default 10,000]
    random_state : None | int
        Random state to use for sampling. Set to None if you want full random
        with each iteration.

    Returns
    -------
    pd.DataFrame
        Rows are the number of UMI sorted cells. Columns are Bulk RNASeq
        samples. Values are Spearman r coefficients.

    """

    res = []
    loc = interval
    while loc < stop:
        idx = umi.sample(n=loc, random_state=random_state).index
        dat = filter_gene_counts_by_barcode(idx, raw).sum(axis=1)
        corrs = []
        for col in bulk_dat.columns:
            corrs.append(spearmanr(bulk_dat[col], dat).correlation)

        res.append([loc, *corrs])
        loc += interval

    col_names = ['Cell Number']
    col_names.extend(bulk_dat.columns)

    df = pd.DataFrame(res, columns=col_names)

    return df.set_index('Cell Number')
