import numpy as np
import scipy.io 
import scrublet as scr

cell_ids = snakemake.input.cell_ids # '../output/cellranger-wf/testis1/outs/raw_gene_bc_matrices/cell_ids.tsv'
mtx = snakemake.input.mtx
cell_calls = snakemake.input.cell_calls

def main():
    cells = get_cell_calls()
    counts_matrix = scipy.io.mmread(mtx)[:, cells].T.tocsc()

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)


def get_cell_calls():
    # read cell_id file
    with open(cell_ids) as fh:
        ids = fh.read().strip().split('\n')

    # read cell_calls file
    pd.read_feather(cell_calls)

    # Create a bool array of cells with 
    flag_cells = []

    return flag_cells


if __name__ == '__main__':
    main()