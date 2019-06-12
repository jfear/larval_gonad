import numpy as np
import scipy.io
import pandas as pd

scrublet = snakemake.input.scrublet
cell_calls = snakemake.input.cell_calls
cell_ids = snakemake.input.cell_ids
genes = snakemake.input.genes
mtx = snakemake.input.mtx
fbgn2symbol_ = snakemake.input.fbgn2symbol

barcodes = snakemake.output.barcodes
genes_out = snakemake.output.genes
mtx_out = snakemake.output.mtx


def main():
    # Remove bad cells from cell ids.
    idx = get_cell_id_index()
    ids = _read_text(cell_ids)[idx]
    with open(barcodes, 'w') as fout:
        fout.write('\n'.join(ids))

    # Remove bad cells from matrix
    mm = scipy.io.mmread(mtx).tocsc()[:, idx].tocoo()
    scipy.io.mmwrite(mtx_out, mm)

    # Add gene symbol
    fbgn2symbol = (pd.read_feather(fbgn2symbol_)
        .set_index('FBgn')
    )

    (pd.read_csv(genes, sep='\t', header=None, index_col=0)
        .rename_axis('FBgn')
        .join(fbgn2symbol, how='left')
        .gene_symbol
        .reset_index()
        .to_csv(genes_out, sep='\t', index=None, header=None)
    )


def _read_text(fname):
    with open(fname) as fh:
        return np.array(fh.read().strip().split('\n'))


def get_intersection_and_dedub():
    df_cell_calls = (
        pd.read_feather(cell_calls)
        .set_index('cell_id')
        # calculate the intersection of cellranger3-wf and droputils
        .loc[:, ['cellranger3-wf', 'droputils']]
        .sum(axis=1) == 2
    )

    dublet_ids = _read_text(scrublet)
    df_cell_calls[df_cell_calls.index.isin(dublet_ids)] = False
    good_cell_ids = df_cell_calls[df_cell_calls].index
    return good_cell_ids


def get_cell_id_index():
    """Get index number for good cells."""
    good_cell_ids = get_intersection_and_dedub()
    matrix_col_names = _read_text(cell_ids)
    return np.in1d(matrix_col_names, good_cell_ids)


if __name__ == "__main__":
    main()
    