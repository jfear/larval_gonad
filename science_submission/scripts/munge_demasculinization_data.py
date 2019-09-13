"""Demasculinization data munging.

Demasculinization plots are barplots with (chrom x % Genes On). Here I munge
data to be easily ploted. The shelve consists of:

* data: [chrom, Proportion of Genes]
* pvalues: [chrom, qval, x, y] where x and y are locations to plot string
  representation.

Plot Example
-------
>>> ax = sns.barplot("chrom", "Proportion of Genes", data=data)
>>> larval_gonad.plotting.stats.format_pval(pvalues.x, pvalues.y, pvalues.qval, ax)

"""
import pandas as pd

from larval_gonad.io import pickle_load, shelve_dump
from larval_gonad.stats import run_chisq


def main():
    fbgn2chrom = (
        pickle_load(SNAKE["fbgn2chrom"]).squeeze().pipe(lambda x: x[x.isin(SNAKE["chrom_order"])])
    )
    num_genes_per_chrom = fbgn2chrom.value_counts()

    # pull out bulk genes based on flag passed {up, down, ns}
    bias = flag_bulk_deg(SNAKE["deg"]).pipe(lambda x: x[x == SNAKE["flag_deg"]])

    num_biased_genes_per_chrom = pd.concat(
        [bias, fbgn2chrom], sort=True, axis=1, join="inner"
    ).chrom.value_counts()

    prop_biased_genes_on_per_chrom = (
        num_biased_genes_per_chrom.div(num_genes_per_chrom)
        .dropna()
        .rename("Proportion of Genes")
        .reindex(SNAKE["chrom_order"])
    )

    # Run stats to determine enrichment
    ct = build_crosstab(num_genes_per_chrom, num_biased_genes_per_chrom)

    qval = (
        run_chisq(ct)
        .loc[("On", "fdr q-value"), :]
        .rename("qval")
        .reindex(SNAKE["chrom_order"])
        .to_frame()
        .rename_axis("chrom")
    )
    qval["x"] = range(qval.shape[0])
    qval["y"] = prop_biased_genes_on_per_chrom
    qval = qval.join(num_biased_genes_per_chrom.rename("chrom_count")).dropna()

    # Save results
    shelve_dump(
        SNAKE["output_file"],
        data=prop_biased_genes_on_per_chrom.to_frame().rename_axis("chrom").reset_index(),
        pvalues=qval,
    )


def flag_bulk_deg(file_name: str, alpha: float = 0.01) -> pd.Series:
    df = pd.read_csv(file_name, sep="\t", index_col=0).assign(bias="ns")
    df.loc[(df.padj <= alpha) & (df.log2FoldChange > 0), "bias"] = "up"
    df.loc[(df.padj <= alpha) & (df.log2FoldChange < 0), "bias"] = "down"
    return df.bias


def build_crosstab(bg_cnt: pd.Series, on_cnt: pd.Series) -> pd.DataFrame:
    """Create crosstab table from backround.

    Uses a background count and an on count to create a crosstab table.

    """
    return (
        pd.concat(
            [on_cnt.rename("On"), (bg_cnt - on_cnt).dropna().rename("Off")], axis=1, sort=True
        )
        .dropna(axis=0)
        .T
    )


def get_direction_com_comparison(gene_set_name):
    """Determine which direction when filtering DEG.

    I always do male/testis vs female/ovary, so male genes are always up and
    female genes are always down. I want to pick this filter based on the
    wildcards provided.

    Example
    -------
    >>> get_direction_com_comparison('male_soma')
    'up'
    >>> get_direction_com_comparison('female_soma')
    'down'
    >>> get_direction_com_comparison('sex_soma_ns')
    'ns'

    """
    parts = gene_set_name.split("_")
    if "ns" in parts:
        return "ns"

    if ("testis" in parts) or ("male" in parts):
        return "up"

    if ("ovary" in parts) or ("female" in parts):
        return "down"


if __name__ == "__main__":
    SNAKE = dict(
        deg=snakemake.input["deg"],
        fbgn2chrom=snakemake.input["fbgn2chrom"],
        output_file=snakemake.output[0],
        chrom_order=snakemake.params["chrom_order"],
    )

    # Debug Settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
    #     print(os.getcwd())
    # except:
    #     pass
    # from larval_gonad.config import read_config
    # config = read_config("../../config/common.yaml")
    # SNAKE = dict(
    #     deg="../../output/expression-atlas-wf/dmel_gonad_biased_expression.tsv",
    #     fbgn2chrom="../../output/x-to-a-wf/fbgn2chrom.pkl",
    #     output_file="",
    #     chrom_order=config["chrom_order"],
    # )

    SNAKE["flag_deg"] = get_direction_com_comparison(snakemake.wildcards.fbgns)

    main()
