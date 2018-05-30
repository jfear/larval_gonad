"""Identification and selection of cells from 10x output.

The first processing step after running cell ranger is determining what
constitutes a cell and what is background. Cell ranger's algorithm does not
seem to work well in our situation due to differences in RNA content. The
number of UMI seems to be the best criteria to perform initial filtering,
perhaps followed by a gene filtering.

This module provides a set of functions for handling parsing of 10x HDF5
formats and munging data do decide what selection criteria should be.

"""
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import tables
from .config import memory

SOMA = [
    'bnb',
    'CadN',
    'cora',
    'crp',
    'dsx',
    'Egfr',
    'egr',
    'ems',
    'Fas3',
    'fax',
    'foxo',
    'fru',
    'gbb',
    'ImpL2',
    'Lar',
    'nord',
    'Nrt',
    'oys',
    'Sox100B',
    'spi',
    'spict',
    'tj',
    'tkv',
    'vkg',
]

EARLY_GERM = [
    'bam',
    'bgcn',
    'hts',
    'Marf'
    'mle',
    'Phf7',
    'Rbp9',
    'tej',
    'tut',
    'vas',
]

LATE_GERM = [
    'aly',
    'c-cup',
    'CG3927',
    'd-cup',
    'fzo',
    'mia',
    'p-cup',
    'r-cup',
    'soti',
    'sowi',
    'sunz',
    'wa-cup',
]

NUCS = ['A', 'C', 'G', 'T']
NUCS_INVERSE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

CellRangerCounts = namedtuple('CellRangerCounts',
                              ['matrix', 'gene_ids', 'barcodes'])


def compress_seq(s: str):
    """ Pack a DNA sequence (no Ns!) into a 2-bit format, in a 64-bit unit

    Most significant bit is set if there was an error

    Based on code from: https://github.com/10XGenomics/cellranger

    cellranger/lib/python/cellranger/utils.py

    """
    bits = 64
    assert len(s) <= (bits/2 - 1)
    result = 0
    for nuc in s:
        if nuc not in NUCS_INVERSE:
            return 1 << (bits - 1)
        result = result << 2
        result = result | NUCS_INVERSE[nuc]
    return result


def decompress_seq(x: int, length=16):
    """ Unpack a DNA sequence from a 2-bit format

    Based on code from: https://github.com/10XGenomics/cellranger

    cellranger/lib/python/cellranger/utils.py

    Parameters
    ----------
    x : int
        Number sequence to be decoded.
    length : int
        Length of the barcode. This can be found in the molecular info hdf5
        file from 10x genome.
        molInfo.get_node_attr('/metrics', 'chemistry_barcode_read_length')

    """
    bits = 64
    x = np.uint64(x)
    assert length <= (bits/2 - 1)
    if x & (1 << (bits-1)):
        return 'N' * length
    result = bytearray(length)
    for i in range(length):
        result[(length-1)-i] = bytearray(NUCS[x & np.uint64(0b11)].encode())[0]
        x = x >> np.uint64(2)
    return result.decode()


def two_bit_mapper(iterable):
    """Return a dictionary mapping 2bit encoded Seqs.

    Parameters
    ----------
    iterable : np.array | pd.Series
        Array of 2bit encoded sequences.

    Returns
    -------
    dict : Mapper from encoded to decoded

    """
    return {k: decompress_seq(k) for k in iterable.unique()}


def decode_cell_names(iterable):
    """Use two_bit_mapper to decode cell names.

    iterable : np.array
        An array of twobit encoded cell names.

    """
    mapper = two_bit_mapper(np.unique(iterable))
    return [mapper[x] for x in iterable]


@memory.cache
def cellranger_umi(fname):
    with tables.open_file(fname, 'r') as f:
        group = f.get_node('/')
        cell_ids = getattr(group, 'barcode').read()
        umi = getattr(group, 'umi').read()
        read_cnts = getattr(group, 'reads').read()

    return pd.DataFrame(dict(
        cell_id=cell_ids,
        umi=umi,
        read_cnt=read_cnts
    ))


@memory.cache
def cellranger_counts(fname, genome='dm6.16'):
    """Import cell ranger counts.

    Cell ranger stores it counts tables in a hdf5 formatted file. This reads
    this file and outputs them as a DataFrame.

    Parameters
    ----------
    fname : str
        Name of hdf5 store.
    genome : str
        Group where data is stored.
    barcodes : list of int
        Encoded barcodes names to filter by

    Returns
    -------
    namedtuple: matrix, gene_ids, barcodes

    """
    with tables.open_file(fname, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, 'genes').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()

    matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
    gene_ids = np.array([x.decode() for x in gene_ids])
    barcodes = np.array([x.decode().replace('-1', '') for x in barcodes])

    return CellRangerCounts(matrix, gene_ids, barcodes)


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


def plot_barcode_rank(background, selected=None, title=None, ax=None):
    """Plot Barcode Rank Plot.

    Example
    -------
    >>> selected = background.query(f'umi_count > {cutoff}')
    >>> barcode_rank_plot_with_cells(background, selected, 'test')

    """
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    ax.plot(range(len(background)), background, label='background')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Unique Cell (log10)')
    ax.set_ylabel('UMI Count (log10)')

    if title is not None:
        ax.set_title(title)

    if selected is not None:
        ax.plot(range(len(selected)), selected, label='cell')

    return ax
