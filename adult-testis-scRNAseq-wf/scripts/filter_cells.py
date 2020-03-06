import os

import numpy as np
import pandas as pd
import scipy.io
from scipy.sparse import coo_matrix


def main():
    good_cells = get_good_cells(snakemake.input.cell_calls, snakemake.input.scrublet)

    tenx_data = load_10x(
        snakemake.input.matrix, snakemake.input.cell_ids, snakemake.input.features, good_cells
    )

    flag_pass_gene_expression_threshold = (
        (tenx_data > 0)
        .sum()
        .map(lambda x: (snakemake.params.low_gene < x) & (x <= snakemake.params.high_gene))
    )

    fbgn2symbol = pd.read_feather(
        snakemake.input.gene_annotation, columns=["FBgn", "gene_symbol"]
    ).set_index("FBgn")

    tenx_data = (
        tenx_data.join(fbgn2symbol, how="left")
        .set_index("gene_symbol", append=True)
        .pipe(lambda df: df.loc[:, flag_pass_gene_expression_threshold])
    )

    # Write output
    ## genes
    tenx_data.index.to_frame().reset_index(drop=True).to_csv(
        snakemake.output.features, sep="\t", index=False, header=False
    )

    ## barcodes
    tenx_data.columns.to_frame().reset_index(drop=True).to_csv(
        snakemake.output.barcodes, sep="\t", index=False, header=False
    )

    ## matrix
    mm = coo_matrix(tenx_data.values)
    scipy.io.mmwrite(snakemake.output.matrix, mm)


def get_good_cells(cell_calls, scrublet):
    non_empty = (
        pd.read_feather(cell_calls, columns=["cell_id", "is_cell"])
        .pipe(lambda df: df[df.is_cell])
        .set_index("cell_id")
        .index
    )

    doublets = pd.read_csv(scrublet, header=None, names=["cell_id"]).set_index("cell_id").index

    return non_empty ^ doublets


def load_10x(mtx, cell_ids, genes, good_cells):
    _barcodes = (
        pd.read_csv(cell_ids, sep="\t", header=None, names=["cell_id"]).set_index("cell_id").index
    )

    _genes = (
        pd.read_csv(genes, sep="\t", header=None, names=["FBgn", "FBgn2"]).set_index("FBgn").index
    )

    flag_good_cells = _barcodes.isin(good_cells)

    # Remove bad cells from matrix
    mm = scipy.io.mmread(mtx).tocsc()[:, flag_good_cells].todense()

    return pd.DataFrame(mm, index=_genes, columns=_barcodes[flag_good_cells])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="adult-testis-scRNAseq-wf",
            input=dict(
                cell_ids="../output/adult-testis-scRNAseq-wf/Ral517/outs/raw_feature_bc_matrix/cell_ids.tsv",
                features="../output/adult-testis-scRNAseq-wf/Ral517/outs/raw_feature_bc_matrix/features.tsv",
                matrix="../output/adult-testis-scRNAseq-wf/Ral517/outs/raw_feature_bc_matrix/matrix.mtx",
                cell_calls="../output/adult-testis-scRNAseq-wf/Ral517_combined_cell_calls.feather",
                scrublet="../output/adult-testis-scRNAseq-wf/Ral517_scrublet_dublets.txt",
                gene_annotation=f"../references/gene_annotation_dmel_r6-24.feather",
            ),
            params=dict(threshold=0.25),
        )

    main()
