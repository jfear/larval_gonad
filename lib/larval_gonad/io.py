"""Collection of io related items."""
import pandas as pd
import scipy.sparse as sp_sparse
import tables


# Import Cell Ranger molecular info
def cellranger_counts(fname, genome='dm6.16'):
    """Import cell ranger counts.

    Cell ranger stores it counts tables in a hdf5 formatted file. This reads
    this file and outputs them as a DataFrame.

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
    return pd.DataFrame(
        data=matrix.todense(),
        index=[x.decode() for x in gene_ids],
        columns=[x.decode() for x in barcodes]
    )
