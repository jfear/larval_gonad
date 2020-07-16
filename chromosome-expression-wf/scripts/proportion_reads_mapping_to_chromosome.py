"""Calculate the proportion of reads mapping to each chromosome element."""
from typing import Generator

import pandas as pd
import joblib
from more_itertools import chunked


CHROMS = ["X", "Y", "2L", "2R", "3L", "3R", "4"]


def main():

    pd.concat(
        [calculate_proportions(raw) for raw in load_gene_level_counts()],
        sort=False,
        ignore_index=True,
    ).to_feather(snakemake.output[0])


def load_gene_level_counts() -> Generator[pd.DataFrame, None, None]:
    """Load gene level UMI counts for each cell.

    Returns:
        pd.DataFrame: FBgn,cell_id,UMI,chrom
    """
    target_fbgns = joblib.load(
        snakemake.input.gene_list
    )  # pylint: disable=unused-variable

    gene_annotation = (
        pd.read_feather(snakemake.input.gene_annotation, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .rename("chrom")
    )

    # For all expressed genes I keep running out memory, so I am going to chunk
    cell_ids = (
        pd.read_feather(snakemake.input.raw).drop("FBgn", axis=1).columns.tolist()
    )
    for chunk in chunked(cell_ids, 1000):
        yield (
            pd.melt(
                pd.read_feather(snakemake.input.raw, columns=["FBgn",] + chunk).query(
                    "FBgn in @target_fbgns"
                ),
                id_vars="FBgn",
                var_name="cell_id",
                value_name="UMI",
            )
            .merge(gene_annotation, on="FBgn")
            .query("chrom in @CHROMS")
        )


def calculate_proportions(raw: pd.DataFrame) -> pd.DataFrame:
    prop_reads_per_chrom = cell_level_prop_reads_mapping_per_chrom(raw)
    scaled_prop_reads_per_chrom = scale_by_number_gene_per_chrom(
        raw, prop_reads_per_chrom
    )

    return pd.concat(
        [prop_reads_per_chrom, scaled_prop_reads_per_chrom], sort=False, axis=1
    ).reset_index()


def cell_level_prop_reads_mapping_per_chrom(raw: pd.DataFrame,) -> pd.Series:
    """Proportion of reads mapping to each chromosome (by cell).

    Args:
        raw (pd.DataFrame): Raw gene level UMI counts.

    Returns:
        pd.Series: [cell_id,chrom],prop_reads
    """
    total_reads_per_cell = raw.groupby("cell_id").UMI.sum()
    reads_per_chrom_per_cell = raw.groupby(["cell_id", "chrom"]).UMI.sum()
    return (reads_per_chrom_per_cell / total_reads_per_cell).rename("prop_reads")


def scale_by_number_gene_per_chrom(
    raw: pd.DataFrame, prop_reads_per_chrom: pd.DataFrame, factor: int = 100
) -> pd.Series:
    """Scale proportions by the number of expressed genes per chromosome.

    Args:
        factor (int, optional): Scaling factor to make the number of genes
        per chromosome close to 1. Defaults to 10.

    Returns:
        pd.Series: [cell_id,chrom],scaled_prop_reads
    """
    scaled_genes_per_chrom_per_cell = (
        raw.groupby(["cell_id", "chrom"]).FBgn.size() / factor
    )

    return (prop_reads_per_chrom / scaled_genes_per_chrom_per_cell).rename(
        "scaled_prop_reads"
    )


if __name__ == "__main__":
    main()
