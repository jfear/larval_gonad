"""Male-Biased Gene Set Enrichment Analysis"""
from more_itertools import chunked
import pandas as pd
import joblib
import gseapy as gp

from larval_gonad.normalization import tpm

THREADS = snakemake.threads


def main():
    expressed_fbgns = joblib.load(snakemake.input.expressed)
    male_biased_fbgs = load_male_biased_genes(expressed_fbgns)

    df = pd.concat(
        [
            # GSEA Results
            pd.concat(
                [
                    run_gsea(cells, male_biased_fbgs)
                    for cells in load_data(expressed_fbgns)
                ],
                sort=False,
            ),
            # Cell ID to Cluster Mapping
            (
                pd.read_feather(
                    snakemake.input.clusters, columns=["cell_id", "cluster"]
                )
                .set_index("cell_id")
                .squeeze()
            ),
        ],
        axis=1,
        sort=False,
        join="inner",
    )

    df.reset_index().to_feather(snakemake.output[0])


def load_male_biased_genes(expressed_genes: list):
    """Load the testis biased genes that are expressed in the single cell data."""
    return set(joblib.load(snakemake.input.male_biased)).intersection(
        set(expressed_genes)
    )


def load_data(expressed_genes: list):
    gene_lengths = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )

    cell_ids = pd.read_feather(snakemake.input.sc).set_index("FBgn").columns.tolist()
    for chunk in chunked(cell_ids, 1000):
        raw = (
            pd.read_feather(snakemake.input.sc, columns=["FBgn"] + chunk)
            .set_index("FBgn")
            .reindex(expressed_genes)
        )
        yield tpm(raw, gene_lengths.reindex(expressed_genes))


def run_gsea(data, gene_set, **kwargs):
    defaults = dict(
        outdir="test/ssgsea_report",
        sample_norm_method="rank",  # choose 'custom' for your own rank list
        permutation_num=0,  # skip permutation procedure, because you don't need it
        no_plot=True,  # skip plotting, because you don't need these figures
        max_size=len(gene_set),
        processes=THREADS,
        format="png",
        seed=9,
    )
    defaults.update(kwargs)

    ss = gp.ssgsea(data=data, gene_sets={"male": gene_set}, **defaults)
    return pd.DataFrame(ss.resultsOnSamples).T.squeeze()


if __name__ == "__main__":
    main()
