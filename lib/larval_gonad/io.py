"""Collection of io related items."""
from pathlib import Path
from joblib import Memory
from yaml import load

import pandas as pd
import scipy.sparse as sp_sparse
import tables

# Create Cache
cache_dir = Path(Path(__file__).parent, '../../output/cache').resolve()
cache_dir.mkdir(exist_ok=True)
memory = Memory(cachedir=cache_dir, verbose=0)

# Get config
cname = Path(Path(__file__).parent, '../../config/common.yml').resolve()
with cname.open() as fh:
    config = load(fh)


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
    return pd.DataFrame(data=matrix, index=gene_ids, columns=barcodes)
