"""Chromosome enrichment tests for different gene subsets."""

import joblib
import pandas as pd

from larval_gonad.stats import run_chisq

CHROMS = snakemake.params.chrom_order


def main():
    all_fbgn2chrom = load_fbgn2chrom()

    expressed = joblib.load(snakemake.input.expressed_fbgns) # pylint: disable=unused-variable
    expressed_fbgn2chrom = all_fbgn2chrom.query("FBgn in @expressed")

    df = pd.concat(
        [
            run_stats(snakemake.input.expressed_fbgns, "All Expressed", all_fbgn2chrom),
            run_stats(snakemake.input.tau_fbgns, "Tau", expressed_fbgn2chrom),
            run_stats(snakemake.input.tsps_fbgns, "TSPS", expressed_fbgn2chrom),
            run_stats(
                snakemake.input.widely_expressed_fbgns,
                "Widely Expressed",
                expressed_fbgn2chrom,
            ),
        ],
        sort=False,
    ).reindex(columns=CHROMS)

    df.columns.name = ""
    df.to_csv(snakemake.output[0], sep="\t")


def load_fbgn2chrom():
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .rename(columns={"FB_chrom": "chrom"})
        .query("chrom in @CHROMS")
    )


def run_stats(file_name: str, name: str, background: pd.DataFrame) -> pd.DataFrame:
    df = background.copy()
    df["in_set"] = False

    # Flag genes in set
    gene_set = joblib.load(file_name)
    df.loc[df.FBgn.isin(gene_set), "in_set"] = True

    return (
        run_chisq(pd.crosstab(df.in_set, df.chrom))
        .loc[(True, slice(None)), :]
        .reset_index()
        .assign(in_set=name)
        .set_index(["in_set", "type"])
    )


if __name__ == "__main__":
    main()
