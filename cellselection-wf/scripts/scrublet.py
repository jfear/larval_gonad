import numpy as np
import scipy.io 
import pandas as pd

import scrublet as scr

cell_ids = snakemake.input.cell_ids 
mtx = snakemake.input.mtx
cell_calls = snakemake.input.cell_calls

DUBS = snakemake.output.doublets
HIST = snakemake.output.histogram
UMAP = snakemake.output.umap
TSNE = snakemake.output.tsne

THRESHOLD = snakemake.params.threshold

def main():
    # read in ids of all columns from the raw count matrix
    with open(cell_ids) as fh:
        cells = np.array(fh.read().strip().split('\n'))

    # Use cellranger and dropletutils to pull out non-empty cells.
    calls = pd.read_feather(cell_calls)
    filtered_cells = calls[calls.sum(axis=1) >= 3].cell_id.values

    # Create a bool array of non-empty cells
    cell_index = np.in1d(cells, filtered_cells)

    # Import raw counts matrix removing empty cells
    counts_matrix = scipy.io.mmread(mtx)[:, cells].tocsc()[:, cells].T

    # Run scrublet
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    flag_doublet = scrub.call_doublets(threshold=THRESHOLD)

    # Save out cell ids that are doublets
    dubs = cells[cell_index][flag_doublet]
    with open(DUBS, 'w') as fout:
        fout.write('\n'.join(dubs))

    # Plot distribution to check threshold
    scrub.plot_histogram()
    fig = plt.gcf()
    fig.savefig(HIST)

    # Plot tSNE embedding
    scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, perplexity=30))
    scrub.plot_embedding('tSNE', order_points=True)
    fig = plt.gcf()
    fig.savefig(TSNE)
    
    # Plot UMAP embedding
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)
    fig = plt.gcf()
    fig.savefig(UMAP)


if __name__ == '__main__':
    main()