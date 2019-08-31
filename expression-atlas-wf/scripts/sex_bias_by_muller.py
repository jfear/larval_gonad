"""Summarize sex biased expression by Muller element.

1. Read differential sex-biased expression results.
2. Read in Muller element mapping.
3. Tabulate the number of genes Up/Down/NS per Muller element.
5. Calculate proportions of Up/Down/NS per Muller element for plotting.
4. Calculate enrichment using chi-sq with post-hoc test.

"""
import os
import pandas as pd

from larval_gonad.io import pickle_load, shelve_dump
from larval_gonad.stats import run_chisq


def main():
    deg = read_deg(snakemake.input.deg, float(snakemake.params.alpha))
    bias_counts_by_muller = count_by_muller(snakemake.input.muller, deg)
    muller_proportions = calc_proprotions(bias_counts_by_muller)
    enrichment = enrichment_analysis(bias_counts_by_muller)

    # Pull out q-values with locations for easy plotting.
    male_qval = pull_out_qval(enrichment, muller_proportions, "male")
    female_qval = pull_out_qval(enrichment, muller_proportions, "female")

    shelve_dump(
        snakemake.output[0], data=muller_proportions, male_qval=male_qval, female_qval=female_qval
    )


def read_deg(file_deg: str, alpha: float = 0.01) -> pd.DataFrame:
    """Read deg and make bias flags.

    Make a flag for male, female, and NS.

    """
    deg = pd.read_csv(file_deg, sep="\t").set_index("YOgn").fillna(1).assign(bias="NS")
    sig_mask = deg.padj <= alpha
    deg.loc[sig_mask & (deg.log2FoldChange > 0), "bias"] = "male"
    deg.loc[sig_mask & (deg.log2FoldChange < 0), "bias"] = "female"

    return deg


def count_by_muller(file_muller: str, deg: pd.DataFrame) -> pd.DataFrame:
    # Read muller element
    yo2muller = pd.Series(pickle_load(file_muller), name="muller").rename_axis("YOgn")

    # Merge and count Male/Female/NS by muller
    return (
        deg.join(yo2muller)
        .groupby("muller")
        .bias.value_counts()
        .unstack()
        .fillna(0)
        .loc[["A", "B", "C", "D", "E", "F"], ["male", "female", "NS"]]
    )


def calc_proprotions(counts_by_muller: pd.DataFrame) -> pd.DataFrame:
    col_order = ["male", "NS", "female"]
    return counts_by_muller.pipe(lambda x: x.div(x.sum(axis=1), axis=0)).loc[:, col_order]


def enrichment_analysis(counts_by_muller: pd.DataFrame) -> pd.DataFrame:
    return (
        run_chisq(counts_by_muller.T)
        .loc[(["male", "female"], "fdr q-value"), :]
        .droplevel("type")
        .T
    )


def pull_out_qval(enrichment: pd.DataFrame, proportions: pd.DataFrame, sex: str):
    return (
        enrichment[sex]
        .rename("q-value")
        .to_frame()
        .assign(x=lambda x: range(x.shape[0]))
        .join(proportions[sex].rename("y"))
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input={
                "muller": "../output/expression-atlas-wf/YOgn_to_muller/dmel.pkl",
                "deg": "../output/expression-atlas-wf/sex_biased_expression/orgR_AC.tsv",
            },
            params=dict(alpha=0.01),
        )

    main()
