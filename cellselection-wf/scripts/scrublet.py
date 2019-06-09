import numpy as np
import scipy.io 
import scrublet as scr

barcodes = snakemake.input.barcodes
mtx = snakemake.input.mtx
cell_calls = snakemake.input.cell_calls

def main():
    cells = get_cell_calls()
    counts_matrix = scipy.io.mmread(mtx)[:, cells].T.tocsc()

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)


def get_cell_calls():

    # read barcodes file

    # read cell_calls file

    # Create a bool array of cells with 
    flag_cells = []

    return flag_cells



if __name__ == '__main__':
    main()