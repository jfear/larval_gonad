"""Merge sex biased expression data set for D. melanogaster.

Takes differential expression results for w1118 and OreR for all tissues and
merges it into a single table. This will allow easier figure production and
summary statistics.
"""
from pathlib import Path

import joblib
import pandas as pd


def main():
    (
        pd.concat([read_table(file_name) for file_name in snakemake.input.deg])
        .pipe(map_to_fbgn)
        .pipe(keep_genes_with_fbgn)
        .pipe(annotate_bias)
        .pipe(rename_orgR_to_OreR)
        .rename(columns={"YOgn": "FBgn"})
        .reset_index(drop=True)
        .to_feather(snakemake.output[0])
    )


def read_table(file_name: str) -> pd.DataFrame:
    """Read in the DESeq2 results table and add species/tissue columns."""
    sample_name = Path(file_name).stem
    species, tissue = sample_name.split("_")
    return (
        pd.read_csv(file_name, sep="\t").assign(species=species).assign(tissue=tissue)
    )


def map_to_fbgn(df: pd.DataFrame) -> pd.DataFrame:
    """Convert YOgn to FBgn where possible."""
    annot = joblib.load(snakemake.input.yogn2fbgn)
    return df.set_index("YOgn").rename(index=annot).reset_index()


def keep_genes_with_fbgn(df: pd.DataFrame) -> pd.DataFrame:
    """Remove genes that do not have an FBgn.

    Haiwang re-annotated the genome, adding new genes that do not have a
    FlyBase annotation. These genes are not of particular importance and
    would just add complexity, so I am removing them.
    """
    return df[df.YOgn.str.startswith("FBgn")]


def annotate_bias(df: pd.DataFrame, alpha: float = 0.01) -> pd.DataFrame:
    """Adds a column annotating the direction of bias (Female, NS, Male)."""
    significant = df.padj <= alpha
    male_bias = df.log2FoldChange > 0
    female_bias = df.log2FoldChange < 0

    df["bias"] = "NS"
    df.loc[significant & male_bias, "bias"] = "male"
    df.loc[significant & female_bias, "bias"] = "female"

    return df


def rename_orgR_to_OreR(df: pd.DataFrame) -> pd.DataFrame:
    df.species = df.species.replace({"orgR": "OreR"})
    return df


if __name__ == "__main__":
    main()
