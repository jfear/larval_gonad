import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io

import scrublet as scr


def main():
    # read in ids of all columns from the raw count matrix
    with open(snakemake.input.cell_ids) as fh:
        cells = np.array(fh.read().strip().split("\n"))

    # Use cellranger and dropletutils to pull out non-empty cells.
    calls = pd.read_feather(snakemake.input.cell_calls, columns=["cell_id", "is_cell"])
    filtered_cells = calls[calls.is_cell].cell_id.values

    # Create a bool array of non-empty cells
    cell_index = np.in1d(cells, filtered_cells)

    # Import raw counts matrix removing empty cells
    counts_matrix = scipy.io.mmread(snakemake.input.matrix).tocsc()[:, cell_index].T

    # Run scrublet
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    _, _ = scrub.scrub_doublets(
        min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
    )
    flag_doublet = scrub.call_doublets(threshold=snakemake.params.threshold)

    # Save out cell ids that are doublets
    dubs = cells[cell_index][flag_doublet]
    with open(snakemake.output.doublets, "w") as fout:
        fout.write("\n".join(dubs))

    # Plot distribution to check threshold
    scrub.plot_histogram()
    fig = plt.gcf()
    fig.savefig(snakemake.output.histogram, dpi=300)

    # Plot tSNE embedding
    scrub.set_embedding("tSNE", scr.get_tsne(scrub.manifold_obs_, perplexity=30))
    scrub.plot_embedding("tSNE", order_points=True)
    fig = plt.gcf()
    fig.savefig(snakemake.output.tsne, dpi=300)

    # Plot UMAP embedding
    scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding("UMAP", order_points=True)
    fig = plt.gcf()
    fig.savefig(snakemake.output.umap, dpi=300)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", None):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="adult-testis-scRNAseq-wf",
            input=dict(
                cell_ids="../output/adult-testis-scRNAseq-wf/Ral517/outs/raw_feature_bc_matrix/cell_ids.tsv",
                matrix="../output/adult-testis-scRNAseq-wf/Ral517/outs/raw_feature_bc_matrix/matrix.mtx",
                cell_calls="../output/adult-testis-scRNAseq-wf/Ral517_combined_cell_calls.feather",
            ),
            params=dict(threshold=0.25)
        )

    main()
