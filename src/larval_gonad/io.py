"""Collection of io related items."""
import os
import shelve
from collections import namedtuple
import pickle

import numpy as np
import pandas as pd
import re
import scipy.sparse as sp_sparse
import tables

NUCS = ["A", "C", "G", "T"]
NUCS_INVERSE = {"A": 0, "C": 1, "G": 2, "T": 3}


CellRangerCounts = namedtuple("CellRangerCounts", ["matrix", "gene_ids", "barcodes"])


def safe_gene_name(symbol):
    """Normalize gene symbols for use as file names."""
    return (
        symbol.replace("(", "")
        .replace(")", "")
        .replace(":", "")
        .replace("&", "")
        .replace("|", "")
        .replace(".", "")
        .replace(" ", "")
    )


def pickle_dump(obj: object, file_name: str):
    with open(file_name, "wb") as handler:
        pickle.dump(obj, handler)


def pickle_load(file_name: str):
    with open(file_name, "rb") as handler:
        return pickle.load(handler)


def shelve_dump(file_name: str, **kwargs):
    """Save a set of objects to a shelve.

    Parameters
    ----------
    file_name: The name of one of the files from a shelve.
    **kwargs: Any number of key word arguments to save to the shelve.

    """
    with shelve.open(os.path.splitext(file_name)[0]) as db:
        for k, v in kwargs.items():
            db[k] = v


def shelve_load(file_name: str, *args):
    """Load a set of objects from a shelve.

    Parameters
    ----------
    file_name: The name of one of the files from a shelve.
    *args: depreciated, does nothing.

    """
    res = {}
    with shelve.open(os.path.splitext(file_name)[0]) as db:
        for k, v in db.items():
            res[k] = v
    return res


def compress_seq(s: str):
    """ Pack a DNA sequence (no Ns!) into a 2-bit format, in a 64-bit uint

    Most significant bit is set if there was an error

    Based on code from: https://github.com/10XGenomics/cellranger

    cellranger/lib/python/cellranger/utils.py

    """
    bits = 64
    assert len(s) <= (bits / 2 - 1)
    result = 0
    for nuc in s:
        if nuc not in NUCS_INVERSE:
            return 1 << (bits - 1)
        result = result << 2
        result = result | NUCS_INVERSE[nuc]
    return result


def decompress_seq(x: int, length=16):
    """ Un-pack a DNA sequence from a 2-bit format

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
    assert length <= (bits / 2 - 1)
    if x & (1 << (bits - 1)):
        return "N" * length
    result = bytearray(length)
    for i in range(length):
        result[(length - 1) - i] = bytearray(NUCS[x & np.uint64(0b11)].encode())[0]
        x = x >> np.uint64(2)
    return result.decode()


def two_bit_mapper(iterable):
    """Return a dictionary mapping 2bit encoded Seqs.

    Parameters
    ----------
    iterable : list-like
        Unique list of 2bit encoded sequences.

    Returns
    -------
    dict : Mapper from encoded to decoded

    """
    return {k: decompress_seq(k) for k in iterable}


def decode_cell_names(iterable):
    """Use two_bit_mapper to decode cell names.

    iterable : np.array
        An array of twobit encoded cell names.

    """
    mapper = two_bit_mapper(np.unique(iterable))
    return [mapper[x] for x in iterable]


def cellranger_umi(fname):
    with tables.open_file(fname, "r") as f:
        group = f.get_node("/")
        barcodes = getattr(group, "barcodes").read()
        barcode_idx = getattr(group, "barcode_idx").read()
        umi = getattr(group, "umi").read()
        read_cnts = getattr(group, "count").read()

    cell_ids = np.char.decode(barcodes[barcode_idx])

    return pd.DataFrame(dict(cell_id=cell_ids, umi=umi, read_cnt=read_cnts))


def cellranger_counts(fname, genome="matrix"):
    """Import cell ranger counts.

    Cell ranger stores it counts tables in a hdf5 formatted file. This reads
    this file and outputs them as a named tuple.

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
    with tables.open_file(fname, "r") as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, "features/id").read()
        barcodes = getattr(group, "barcodes").read()
        data = getattr(group, "data").read()
        indices = getattr(group, "indices").read()
        indptr = getattr(group, "indptr").read()
        shape = getattr(group, "shape").read()

    matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
    gene_ids = np.array([x.decode() for x in gene_ids])
    barcodes = np.array([x.decode().replace("-1", "") for x in barcodes])

    return CellRangerCounts(matrix, gene_ids, barcodes)


def feather_to_cluster_rep_matrix(fname):
    """Helper function to building a cluster rep matrix from a feather"""
    return (
        pd.read_feather(fname)
        .set_index(["FBgn", "cluster", "rep"])
        .iloc[:, 0]
        .unstack(level=[1, 2])
    )


def feather_to_cluster_matrix(fname):
    """Helper function to building a cluster matrix from a feather"""
    return (
        pd.read_feather(fname)
        .set_index(["FBgn", "cluster"])
        .iloc[:, 0]
        .unstack()
    )


def melt_cluster_rep_matrix(df, name="count"):
    """Helper function to melt a cluster rep matrix for saving as feather"""
    return df.T.reset_index().melt(id_vars=["cluster", "rep"], var_name="FBgn", value_name=name)
    

def melt_cluster_matrix(df, name="count"):
    """Helper function to melt a cluster matrix for saving as feather"""
    return df.T.reset_index().melt(id_vars="cluster", var_name="FBgn", value_name=name)



class GffRow(object):
    def __init__(self, row):
        self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes = row.strip().split(
            "\t"
        )
        self.is_gene = self.type == "gene"
        self.parsed_attributes = self.parse_attributes()

    def parse_attributes(self):
        parsed_attributes = {}
        for attr in self.attributes.split(";"):
            mm = re.search('(?P<key>.*?)\s+"(?P<value>.*?)"', attr)
            if mm:
                parsed_attributes[mm.group("key").strip()] = mm.group("value").strip()
        return parsed_attributes

    def __getitem__(self, key):
        return self.parsed_attributes[key]


if __name__ == "__main__":
    main()
