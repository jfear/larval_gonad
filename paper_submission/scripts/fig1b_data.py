"""Combine Adult and Bulk Sex-Biased Expression Data."""

import pandas as pd

CHROMS = [
    "X",
    "Y",
    "2L",
    "2R",
    "3L",
    "3R",
    "4",
]


def main():
    df = (
        load_deg_data(snakemake.input.adult, snakemake.input.larval)
        .pipe(add_chromosome_annotation, snakemake.input.gene_annot)
        .query("chrom in @CHROMS")
        .reset_index(drop=True)
    )

    df.to_feather(snakemake.output[0])


def load_deg_data(adult: str, larval: str) -> pd.DataFrame:
    return pd.concat(
        [pd.read_feather(adult), read_larval_data(larval)],
        sort=False,
        ignore_index=True,
    )


def add_chromosome_annotation(df: pd.DataFrame, gene_annot: str) -> pd.DataFrame:
    """Merge chromosome annotations"""
    annot = pd.read_feather(gene_annot, columns=["FBgn", "FB_chrom"]).rename(
        columns={"FB_chrom": "chrom"}
    )
    return df.merge(annot)


def assign_bias(df: pd.DataFrame, alpha: float = 0.01) -> pd.DataFrame:
    sig = df.padj <= alpha
    male = df.log2FoldChange > 0
    female = df.log2FoldChange < 0

    df["bias"] = "NS"
    df.loc[sig & male, "bias"] = "male"
    df.loc[sig & female, "bias"] = "female"
    return df


def read_larval_data(tsv: str) -> pd.DataFrame:
    """Wrangle larval data to match bulk."""
    return (
        pd.read_csv(tsv, sep="\t")
        .assign(species="w1118", tissue="L3_GO")
        .pipe(assign_bias)
        .drop("stat", axis=1)
    )


if __name__ == "__main__":
    main()
