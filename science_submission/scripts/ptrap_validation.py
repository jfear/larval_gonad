"""Cell type validation code.

This module houses the cell type validation code.
"""
from functools import partial
import numpy as np
import pandas as pd

from larval_gonad.validation import GeneValidator


def main():
    fbgn2symbol = (
        pd.read_feather(snakemake.input["gene_annot"], columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .squeeze()
    )
    ptraps = get_ptraps()
    biomarkers = get_biomarkers()
    zscores = get_zscores()

    protein = True
    res = []
    for (fbgn, ptrap_id), ptrap_gene in ptraps.items():
        gene_score = GeneValidator(
            fbgn,
            ptrap_gene,
            set(),
            biomarkers.get(fbgn, set()),
            zscores.get(fbgn, set()),
            flag_protein=protein,
        )
        res.append(
            [
                fbgn,
                ptrap_id,
                gene_score.lit_gene,
                gene_score.biomarker,
                gene_score.zscore,
                gene_score.score,
            ]
        )
    df = (
        pd.DataFrame(
            res, columns=["FBgn", "ptrap_id", "ptrap_cell_types", "Biomarker", "Zscore", "Score"]
        )
        .fillna(0)
        .set_index("FBgn")
        .join(fbgn2symbol)
        .set_index("gene_symbol", append=True)
    )
    df.to_csv(snakemake.output["table"], sep="\t")
    df.Score.value_counts().rename("Counts").to_frame().to_csv(
        snakemake.output["summary"], sep="\t"
    )


def cell_type_above_upper_quantile(x, n=3, name="zscore"):
    x["flag_upper"] = x[name] > np.quantile(x[name], 0.75)
    flag_cluster = x.groupby(["cluster"]).flag_upper.sum() == n
    return set(flag_cluster[flag_cluster].index.tolist())


def get_ptraps():
    df = (
        pd.read_csv(snakemake.input.ptraps, sep="\t", na_values='-')
        .fillna(0)
        .set_index(["FBgn", "gene_symbol", "ptrap_id"])
    )

    return (
        pd.DataFrame(
            [(*idx, set(dd[dd == dd.max()].index.values.tolist())) for idx, dd in df.iterrows()],
            columns=["FBgn", "gene_symbol", "ptrap_id", "cell_types"],
        )
        .set_index(["FBgn", "ptrap_id"])
        .cell_types.squeeze()
    )


def get_biomarkers():
    return (
        pd.read_feather(snakemake.input.biomarkers)
        .query("p_val_adj <= 0.01")
        .assign(cluster=lambda x: x.cluster.map(snakemake.params["cluster_annot"]))
        .groupby("FBgn")
        .apply(lambda x: set(x.cluster))
        .rename("sc")
    )


def get_zscores():
    return (
        pd.read_feather(snakemake.input.zscores)
        .groupby("FBgn")
        .apply(partial(cell_type_above_upper_quantile, n=3, name="zscore"))
        .rename("zscore")
    )


if __name__ == "__main__":
    main()
