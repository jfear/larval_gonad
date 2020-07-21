"""Run pairwise comparison of enrichment scores."""
import pandas as pd
import seaborn as sns

from larval_gonad.stats import PairwisePermutationTest

THREADS = snakemake.threads

def main():
    df = pd.read_feather(snakemake.input[0])
    pw = PairwisePermutationTest("cluster", "male", data=df, threads=THREADS)
    pw.fit()

    (
        pw.results
        .assign(
            name1=lambda x: pd.Categorical(x.name1, categories=snakemake.params.order)
        )
        .assign(
            name2=lambda x: pd.Categorical(x.name2, categories=snakemake.params.order)
        )
        .sort_values(["name1", "name2"])
        .set_index(["name1", "name2"])
        .to_csv(snakemake.output[0], sep="\t")
    )



if __name__ == "__main__":
    main()
