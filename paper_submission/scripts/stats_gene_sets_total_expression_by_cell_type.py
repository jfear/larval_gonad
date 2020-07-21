"""Plot total expression (TPM) for different gene sets."""
import joblib
import pandas as pd

from larval_gonad.stats import PairwisePermutationTest

THREADS = snakemake.threads


def main():
    tpm = pd.read_feather(snakemake.input.tpm)
    df = pd.concat(
        [
            create_gene_subset(tpm, snakemake.input.expressed_fbgns, "All Expressed"),
            create_gene_subset(
                tpm, snakemake.input.widely_expressed_fbgns, "Widely Expressed"
            ),
            create_gene_subset(tpm, snakemake.input.tau_fbgns, "Tau"),
            create_gene_subset(tpm, snakemake.input.tsps_fbgns, "TSPS"),
        ],
        ignore_index=True,
    )

    df.groupby("gene_set").apply(run_stats).to_csv(snakemake.output[0], sep="\t")


def create_gene_subset(df: pd.DataFrame, gene_set: str, name: str):
    fbgns = joblib.load(gene_set)  # pylint: disable=unused-variable
    return df.query("FBgn in @fbgns").assign(gene_set=name)


def run_stats(df: pd.DataFrame) -> pd.DataFrame:
    pw = PairwisePermutationTest("cluster", "TPM", data=df, threads=THREADS).fit()
    return (
        pw.results.assign(
            name1=lambda x: pd.Categorical(x.name1, categories=snakemake.params.order)
        )
        .assign(
            name2=lambda x: pd.Categorical(x.name2, categories=snakemake.params.order)
        )
        .sort_values(["name1", "name2"])
        .set_index(["name1", "name2"])
    )


if __name__ == "__main__":
    main()
