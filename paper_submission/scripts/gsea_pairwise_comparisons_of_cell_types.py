"""Run pairwise comparison of enrichment scores."""
import pandas as pd
import seaborn as sns

from larval_gonad.stats import PairwisePermutationTest

THREADS = snakemake.threads

def main():
    df = pd.read_feather(snakemake.input[0])
    pw = PairwisePermutationTest("cluster", "male", data=df, threads=THREADS)
    pw.fit()
    pw.results.to_csv(snakemake.output[0], sep="\t", index=False)



if __name__ == "__main__":
    main()
