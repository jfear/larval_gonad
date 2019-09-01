"""Calculate tsps for each YOgn by sex.

I am calculating tsps for each YOgn by sex separately. A score of ≥1 is highly
tissue specific and a score of 0 is a housekeeping gene. NaN values are genes 
that were not expressed.

"""
import os
import numpy as np
import pandas as pd
from scipy.stats import entropy

from larval_gonad.normalization import tpm


def main():
    # Read in counts and aggregate by gene*tissue*sex
    cnts = pd.read_feather(snakemake.input.counts).set_index("YOgn")
    agg_cnts = aggregate_counts_by_sex_tissue(cnts)

    # TPM normalize counts
    gene_lengths = get_gene_lengths(snakemake.input.metadata)
    norm_cnts = tpm(agg_cnts, gene_lengths).dropna()

    # Calculate tsps scores for males and females
    tsps_scores = pd.concat(
        [
            norm_cnts["m"].apply(tsps, axis=1).rename("male"),
            norm_cnts["f"].apply(tsps, axis=1).rename("female"),
        ],
        axis=1,
        sort=True,
    )

    tsps_scores.reset_index().to_feather(snakemake.output[0])


def aggregate_counts_by_sex_tissue(cnts):
    metadata = cnts.columns.str.extract(
        r"(?P<species>\w+)_(?P<tissue>\w+)_(?P<sex>\w+)_(?P<rep>\w+)"
    )
    metadata.index = cnts.columns
    metadata.index.name = "samplename"

    return (
        cnts.reset_index()
        .melt(id_vars="YOgn", var_name="samplename", value_name="cnt")
        .merge(metadata, left_on="samplename", right_index=True)
        .groupby(["YOgn", "tissue", "sex"])
        .cnt.sum()
        .unstack(level=[-1, -2])
    )


def get_gene_lengths(file_name):
    return (
        pd.read_feather(file_name)
        .set_index("YOgn")
        .assign(gene_length=lambda df: df.end.astype(int) - df.start.astype(int))
        .gene_length
    )


def tsps(x: np.ndarray):
    """Calculate tissue specificity score as defined in:

    > Ravasi, Timothy, Harukazu Suzuki, Carlo Vittorio Cannistraci, Shintaro
    > Katayama, Vladimir B. Bajic, Kai Tan, Altuna Akalin, et al. 2010. “An
    > Atlas of Combinatorial Transcriptional Regulation in Mouse and Man.” Cell
    > 140 (5): 744–52.

    Example
    -------
    >>> tsps(np.array([1, 1, 1]))
    0.0
    >>> tsps(np.array([0, 0, 0]))
    np.nan
    >>> tsps(np.array([1, 0, 0]))
    1.5849625007211563

    """
    if x.sum() == 0:
        return np.nan
    _n = x.shape[0]
    _q = np.array([1 / _n] * _n)
    return entropy(x / x.sum(), _q, 2)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input=dict(
                counts="../output/expression-atlas-wf/aggregated_counts_table/orgR.feather",
                metadata="../output/expression-atlas-wf/YOgn_metadata/dmel.feather",
            ),
        )

    main()
