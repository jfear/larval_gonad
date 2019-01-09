"""Identification and selection of cells from 10x output.

The first processing step after running cell ranger is determining what
constitutes a cell and what is background. Cell ranger's algorithm does not
seem to work well in our situation due to differences in RNA content. The
number of UMI seems to be the best criteria to perform initial filtering,
perhaps followed by a gene filtering.

This module provides a set of functions for handling parsing of 10x HDF5
formats and munging data do decide what selection criteria should be.

"""
import numpy as np
import pandas as pd

from .config import memory
from .plotting import make_ax
from .io import CellRangerCounts, cellranger_counts, cellranger_umi


@memory.cache
def build_umi_gene_count_table(cr_raw, cr_umi):
    cr = cellranger_counts(cr_raw)
    umi = cellranger_umi(cr_umi)

    num_genes_on = calc_num_genes_on(cr)
    umi_cnts = umi.query('read_cnt > 0').groupby('cell_id').size().to_frame()
    umi_cnts.columns = ['umi_cnt']

    return umi_cnts.join(num_genes_on).sort_values('umi_cnt', ascending=False)


def calc_num_genes_on(cr: CellRangerCounts):
    """Number of genes that have >0 reads.

    Series containng the number of genes expressed by cell ID given a
    CellRangerCounts object.

    """
    num_genes_on = np.asarray((cr.matrix > 0).sum(axis=0))[0]
    idx = pd.Index(cr.barcodes, name='cell_id')
    return pd.Series(data=num_genes_on, index=idx, name='gene_cnt')


def cells_with_min_gene_expression(cr: CellRangerCounts, cutoff=200):
    """Get cell ids with more than `cutoff` genes.

    We are tyring to determine what kind of cell filtering should be preformed.
    In downstream analysis we use a gene level filter to remove cells that have
    fewer than `cutoff` expressed genes.  I am trying to determine if the cell
    ranger filter is having any affect on this or not.

    Parameters
    ----------
    cr : CellRangerCounts
        The raw_gene_bc_matrices_h5.h5 from cell ranger.
    cutoff : int
        The minimum number of expressed genes.

    Example
    -------
    >>> cr = cellranger_counts(...)
    >>> cells_on = cells_with_min_gene_expression(cr, cutoff=300)

    """
    num_genes_on = calc_num_genes_on(cr)
    cells_on = num_genes_on > cutoff
    return num_genes_on[cells_on].index.unique().tolist()


def filter_gene_counts_by_barcode(barcodes: np.array, cr: CellRangerCounts):
    """Given the raw gene counts matrix pull out specific cells.

    Example
    -------
    >>> cr = cellranger_counts(...)
    >>> barcodes = [...]
    >>> filter_gene_counts_by_barcode(barcodes, cr)

    """
    mask = np.in1d(cr.barcodes, barcodes)
    bcs = cr.barcodes[mask]
    matrix = cr.matrix[:, mask].todense()

    return pd.DataFrame(
        data=matrix,
        index=pd.Index(cr.gene_ids, name='FBgn'),
        columns=bcs
    )


def get_number_of_expressed_genes(fname):
    """Get number of genes with >0 reads.

    fname : str
        Path to the gene counts table.
    """
    cnts = pd.read_csv(fname, sep='\t')
    return (cnts.sum(axis=1) > 0).sum()


def plot_barcode_rank(umi, selected=None, title=None, **kwargs):
    """Plot Barcode Rank Plot."""
    options = {
        'kind': 'scatter',
        's': 3,
        'logx': True,
        'logy': True
    }
    options.update(kwargs)
    try:
        ax = options.pop('ax')
    except KeyError:
        ax = make_ax()

    dat = umi.to_frame()
    dat.columns = ['umi']

    dat = dat.sort_values('umi', ascending=False)
    dat['cell_num'] = list(range(1, dat.shape[0] + 1))

    dat.plot('cell_num', 'umi', c='lightgrey', ax=ax, **options)

    if selected is not None:
        dat.loc[dat.index.isin(selected), :].plot('cell_num',
                                                  'umi',
                                                  c='g',
                                                  ax=ax,
                                                  **options
                                                  )

    if title is not None:
        ax.set_title(title)
