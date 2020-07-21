"""Run pairwise comparison of enrichment scores."""
import pandas as pd
import seaborn as sns

from larval_gonad.stats import PairwisePermutationTest

THREADS = snakemake.threads


def main():
    df = pd.read_feather(snakemake.input[0])
    pw = PairwisePermutationTest(
        "cluster", "male", data=df, threads=THREADS, order=snakemake.params.order
    ).fit()

    (
        pw.results.sort_values(["name1", "name2"])
        .set_index(["name1", "name2"])
        .to_csv(snakemake.output[0], sep="\t")
    )


if __name__ == "__main__":
    main()
