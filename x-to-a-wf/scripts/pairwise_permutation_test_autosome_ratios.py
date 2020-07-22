import pandas as pd

from larval_gonad.stats import PairwisePermutationTest

THREADS = snakemake.threads


def main():
    df = pd.read_feather(snakemake.input[0])

    results = pd.concat([
        run_stats(df, "x_to_a_ratio"),
        run_stats(df, "fourth_to_a_ratio"),
        run_stats(df, "y_to_a_ratio"),
    ], sort=False)

    results.to_csv(snakemake.output[0], sep="\t")


def run_stats(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Run Pairwise Stats on a given column."""
    pw = PairwisePermutationTest(
        "cluster",
        col,
        df,
        order=snakemake.params.cluster_order,
        threads=THREADS,
    ).fit()

    return pw.results.assign(chrom_ratio=col).set_index(
        ["chrom_ratio", "name1", "name2"]
    )


if __name__ == "__main__":
    main()
