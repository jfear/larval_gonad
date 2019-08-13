"""Calculate autosome ratios for each cell.

This script calculates the {X, 4, and Y} to autosome ratios for each
individual cell. I consider chromosomes 2L, 2R, 3L, and 3R as autosomes.

1. Pull out target FBgns.
2. Sum the number of raw reads for each chromosome.
3. Normalize totals by the number of genes on each chromosome.
4. Take the ratio of X / A, 4 / A, and Y / A

"""
import pandas as pd

from larval_gonad.io import pickle_load


def main(snake):
    annot = gene_annotation_for_target_genes(snake["fbgn2chrom"], snake["target_fbgns"])
    clusters = pd.read_feather(snake["clusters"]).set_index("cell_id")
    num_genes_per_chrom = calculate_number_of_genes_per_chrom(annot, snake["autosomes"])
    agg_counts = aggregate_count_data_to_chrom(snake["raw"], annot, snake["chrom_order"])
    ratios = calculate_ratios(agg_counts, num_genes_per_chrom, snake['autosomes'])

    ratios.join(clusters, how="inner").reset_index().to_feather(snake["output_file"])


def gene_annotation_for_target_genes(fbgn2chrom: str, target_fbgns: str) -> pd.DataFrame:
    """Subset fbg2chrom based on target gene set."""
    return pickle_load(fbgn2chrom).reindex(pickle_load(target_fbgns)).dropna().squeeze()


def calculate_number_of_genes_per_chrom(annot: pd.DataFrame, autosomes: list) -> pd.Series:
    """Count the number of genes on each chromosome and the autosomes together."""
    num_genes_per_chrom = annot.value_counts()
    num_genes_per_chrom["autosome"] = num_genes_per_chrom.loc[autosomes].sum()
    return num_genes_per_chrom


def aggregate_count_data_to_chrom(raw: str, annot: pd.DataFrame, chrom_order: list) -> pd.DataFrame:
    """Sum the number of reads for each chromosome."""
    return (
        pd.read_feather(raw)
        .set_index("FBgn")
        .join(annot, how="inner")
        .groupby("chrom")
        .sum()
        .reindex(chrom_order)
        .fillna(0)
        .T.rename_axis("cell_id")
    )


def calculate_ratios(
    agg_counts: pd.DataFrame, num_genes_per_chrom: pd.Series, autosomes: list
) -> pd.Series:
    """Normalize by gene count and calculate autosome ratios."""
    return (
        agg_counts.assign(autosome=lambda agg_counts: agg_counts[autosomes].sum(axis=1))
        .div(num_genes_per_chrom / 1e3, axis="columns")
        .assign(x_to_a_ratio=lambda agg_counts: agg_counts["X"] / agg_counts.autosome)
        .assign(fourth_to_a_ratio=lambda agg_counts: agg_counts["4"] / agg_counts.autosome)
        .assign(y_to_a_ratio=lambda agg_counts: agg_counts["Y"] / agg_counts.autosome)
        .loc[:, ["x_to_a_ratio", "fourth_to_a_ratio", "y_to_a_ratio"]]
    )


if __name__ == "__main__":
    SNAKE = dict(
        raw=snakemake.input["raw"],
        fbgn2chrom=snakemake.input["fbgn2chrom"],
        clusters=snakemake.input["clusters"],
        target_fbgns=snakemake.input["target_fbgns"],
        autosomes=snakemake.params["autosomes"],
        chrom_order=snakemake.params["chrom_order"],
        output_file=snakemake.output[0],
    )

    # Debug Settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), "x-to-a-wf/scripts"))
    #     print(os.getcwd())
    # except:
    #     pass
    # from larval_gonad.config import read_config
    # config = read_config("../../config/common.yaml")
    # SNAKE = dict(
    #     raw="../../output/cellselection-wf/raw.feather"
    #     fbgn2chrom="../../output/x-to-a-wf/fbgn2chrom.pkl"
    #     clusters="../../output/seurat3-cluster-wf/combined_n3_clusters.feather"
    #     target_fbgns='../../output/cellselection-wf/commonly_expressed_genes.pkl'
    #     snake_autosomes=config["autosomes"]
    #     snake_chrom_order=config["chrom_order"]
    #     snake_output_file=''
    # )

    main(SNAKE)
