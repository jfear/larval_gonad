"""Pairwise comparisons Chromosome Expression By Cell Type."""
import pandas as pd

from larval_gonad.stats import PairwisePermutationTest

THREADS = snakemake.threads


def main():
    df = pd.concat(
        [
            load_data(snakemake.input.all_genes, "All Expressed"),
            load_data(snakemake.input.tau_genes, "Tau"),
            load_data(snakemake.input.tsps_genes, "TSPS"),
            load_data(snakemake.input.widely_expressed_genes, "Widely Expressed"),
        ],
        sort=False,
        ignore_index=True,
    ).merge(pd.read_feather(snakemake.input.clusters), on="cell_id")

    df.groupby("gene_set").apply(run_stats).to_csv(snakemake.output[0], sep="\t")


def load_data(file_name: str, name: str):
    return pd.read_feather(file_name).assign(gene_set=name)


def run_stats(df: pd.DataFrame) -> pd.DataFrame:
    pw = PairwisePermutationTest(
        "chrom",
        "scaled_prop_reads",
        data=df,
        threads=THREADS,
        order=snakemake.params.chrom_order,
    ).fit()
    return pw.results.sort_values(["name1", "name2"]).set_index(["name1", "name2"])


if __name__ == "__main__":
    main()
