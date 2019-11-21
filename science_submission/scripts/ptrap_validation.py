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
    for fbgn, ptrap_gene in ptraps.items():
        gene_score = GeneValidator(
            fbgn,
            ptrap_gene,
            set(),
            biomarkers.get(fbgn, set()),
            zscores.get(fbgn, set()),
            flag_protein=protein,
        )
        res.append(
            [fbgn, gene_score.lit_gene, gene_score.biomarker, gene_score.zscore, gene_score.score]
        )
    df = (
        pd.DataFrame(res, columns=["FBgn", "literature", "Biomarker", "Zscore", "Score"])
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
    return (
        pd.read_csv(snakemake.input.ptraps, sep="\t", na_values='-')
        .fillna(0)
        .groupby("FBgn").mean()
        # Fix column names
        .assign(G=lambda x: x['Spermatogoina'])
        .assign(EPS=lambda x: x['Early 1° Spermatocytes'])
        .assign(MPS=lambda x: x['Mid 1° Spermatocytes'])
        .assign(LPS=lambda x: x['Late 1° Spermatocytes'])
        .assign(C1=lambda x: x['Cyst Cells'])
        .assign(C2=lambda x: x['Cyst Cells'])
        .assign(C3=lambda x: x['Cyst Cells'])
        .assign(C4=lambda x: x['Cyst Cells'])
        .assign(P=lambda x: x['Pigment Cells'])
        .assign(T=lambda x: x['Terminal Epithelium'])
        .loc[:, ["G", "EPS", "MPS", "LPS", "C1", "C2", "C3", "C4", "P", "T"]]
        .reset_index()
        # Figure out what cell types are in upper quantile
        .melt(id_vars="FBgn", var_name="cluster", value_name="score")
        .groupby("FBgn")
        .apply(partial(cell_type_above_upper_quantile, n=1, name="score"))
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
