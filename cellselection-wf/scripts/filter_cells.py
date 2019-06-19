import numpy as np
import scipy.io
from scipy.sparse import coo_matrix
import pandas as pd

SAMPLE = snakemake.wildcards.sample
CELL_IDS = snakemake.input.cell_ids
GENES = snakemake.input.genes
MTX = snakemake.input.mtx
CELL_CALLS = snakemake.input.cell_calls
SCRUBLET = snakemake.input.scrublet
GENE_ANNOTATION = snakemake.input.gene_annotation

BARCODES = snakemake.output.barcodes
GENES_OUT = snakemake.output.genes
MTX_OUT = snakemake.output.mtx

LOW_GENE = snakemake.params.low_gene
HIGH_GENE = snakemake.params.high_gene


def main():
    good_cells = get_good_cells(CELL_CALLS, SCRUBLET)

    tenx_data = load_10x(MTX, CELL_IDS, GENES, good_cells)

    flag_pass_gene_expression_threshold = (
        (tenx_data > 0).sum().map(lambda x: (LOW_GENE < x) & (x <= HIGH_GENE))
    )

    fbgn2symbol = pd.read_feather(
        GENE_ANNOTATION, columns=["FBgn", "gene_symbol"]
    ).set_index("FBgn")

    tenx_data = (
        tenx_data.join(fbgn2symbol, how="left")
        .set_index("gene_symbol", append=True)
        .pipe(lambda df: df.loc[:, flag_pass_gene_expression_threshold])
    )

    # Write output
    ## genes
    tenx_data.index.to_frame().reset_index(drop=True).to_csv(
        GENES_OUT, sep="\t", index=False, header=False
    )

    ## barcodes
    tenx_data.columns.to_frame().reset_index(drop=True).to_csv(
        BARCODES, sep="\t", index=False, header=False
    )

    ## matrix
    mm = coo_matrix(tenx_data.values)
    scipy.io.mmwrite(MTX_OUT, mm)


def get_good_cells(cell_calls, scrublet):
    non_empty = (
        pd.read_feather(cell_calls, columns=["cell_id", "is_cell"])
        .pipe(lambda df: df[df.is_cell])
        .set_index("cell_id")
        .index
    )

    doublets = (
        pd.read_csv(scrublet, header=None, names=["cell_id"]).set_index("cell_id").index
    )

    return non_empty ^ doublets


def load_10x(mtx, cell_ids, genes, good_cells):
    _barcodes = (
        pd.read_csv(cell_ids, sep="\t", header=None, names=["cell_id"])
        .set_index("cell_id")
        .index
    )

    _genes = (
        pd.read_csv(genes, sep="\t", header=None, names=["FBgn", "FBgn2"])
        .set_index("FBgn")
        .index
    )

    flag_good_cells = _barcodes.isin(good_cells)

    # Remove bad cells from matrix
    mm = scipy.io.mmread(mtx).tocsc()[:, flag_good_cells].todense()

    return pd.DataFrame(mm, index=_genes, columns=_barcodes[flag_good_cells])


if __name__ == "__main__":
    main()

