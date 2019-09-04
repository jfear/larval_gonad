import os
import pandas as pd

from larval_gonad.io import shelve_dump
from larval_gonad.stats import run_chisq

MULLER = {"X": "A", "2L": "B", "2R": "C", "3L": "D", "3R": "E", "4": "F"}

def main():
    deg = read_deg(snakemake.input.deg)
    fbgn2muller = (
        pd.read_feather(snakemake.input.annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn").squeeze()
        .map(MULLER).dropna().rename("muller")
    )
    bias_counts_by_muller = count_by_muller(deg.join(fbgn2muller))
    agg_counts_by_muller = aggregate_autosomes(bias_counts_by_muller)
    muller_proportions = calc_proprotions(bias_counts_by_muller)
    enrichment = enrichment_analysis(agg_counts_by_muller)

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
    deg = pd.read_csv(file_deg, sep="\t").set_index("FBgn").fillna(1).assign(bias="NS")
    sig_mask = deg.padj <= alpha
    deg.loc[sig_mask & (deg.log2FoldChange > 0), "bias"] = "male"
    deg.loc[sig_mask & (deg.log2FoldChange < 0), "bias"] = "female"

    return deg


def count_by_muller(deg) -> pd.DataFrame:
    # Merge and count Male/Female/NS by muller
    return (
        deg
        .groupby("muller")
        .bias.value_counts()
        .unstack()
        .fillna(0)
        .loc[["A", "B", "C", "D", "E", "F"], ["male", "female", "NS"]]
    )


def aggregate_autosomes(counts_by_muller):
    """Add autosome counts together."""
    non_autosome = counts_by_muller.reindex(["A", "F"])
    autosomes = (
        counts_by_muller.reindex(["B", "C", "D", "E"]).sum().rename("autosomes").to_frame().T
    )
    return pd.concat([non_autosome, autosomes])


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
    df = (
        enrichment[sex]
        .rename("q-value")
        .to_frame()
        .join(proportions[sex].rename("y"))
        .reindex(["A", "F"])
    )
    df["x"] = 0
    df.loc["F", "x"] = 5
    return df


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="bulk2-rnaseq-wf",
            input=dict(
                deg="../output/bulk2-rnaseq-wf/deg/bulk_testis_vs_ovary.tsv",
                annot="../references/gene_annotation_dmel_r6-26.feather"
            ),
            params=dict(alpha=0.01),
        )

    main()
