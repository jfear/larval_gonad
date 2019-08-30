"""Determine chromosome consensus for scaffolds.

We want to include the out groups D. willistoni, D. virilis, and D.
mojavensus in our analysis. However, genes are assigned to scaffolds and not
Muller elements. We assign chromosome arm based on the consensus from D.
melanogaster.

"""
import numpy as np
import pandas as pd

from larval_gonad.io import pickle_load

DMEL_MULLER = {"X": "A", "2L": "B", "2R": "C", "3L": "D", "3R": "E", "4": "F"}
DPSE_MULLER = {"XL": "A", "4": "B", "3": "C", "XR": "D", "2": "E", "5": "F"}


def main():
    fbgn2muller = (
        pd.read_feather(snakemake.input["gene_annotation"], columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .map(lambda x: DMEL_MULLER.get(x, "unknown_scaffold"))
        .rename("muller")
    )

    primary2secondary = pickle_load(snakemake.input["primary2secondary"])

    orthologs = pd.read_feather(snakemake.input["orthologs"]).set_index("FBgn")

    # update to current reference FBgns
    orthologs.index = orthologs.index.map(primary2secondary)

    # update dmel chromosome based on the annotation
    orthologs = orthologs[~orthologs.index.isnull()].join(fbgn2muller)

    # Sanity check our consensus method with dpse annotation from FlyBase
    dpse_fbgn2annotated = (
        orthologs.dpse.str.extract("([a-zA-Z0-9]+)").squeeze().map(DPSE_MULLER).rename("annotated")
    )

    dpse_consensus = orthologs.groupby("dpse").apply(consensus_muller).rename("muller")
    dpse_fbgn2consensus = orthologs.dpse.map(dpse_consensus).rename("consensus")

    dpse_sanity = pd.concat([dpse_fbgn2annotated, dpse_fbgn2consensus], axis=1).dropna(
        subset=["annotated"]
    )

    # Assure that dpse annotation and consensus are identical
    assert (dpse_sanity.annotated.fillna("NAN") == dpse_sanity.consensus.fillna("NAN")).all()

    # Get muller for other samples
    others = [
        orthologs[x].map(orthologs.groupby(x).apply(consensus_muller).rename("muller"))
        for x in orthologs.columns
        if (x != "dmel") & (x != "dpse") & (x != "muller")
    ]

    df_muller = pd.concat(
        [orthologs.muller.rename("dmel"), dpse_fbgn2consensus.rename("dpse"), *others], axis=1
    )

    df_muller.reset_index().to_feather(snakemake.output[0])


def consensus_muller(x):
    if x.shape[0] > 20:
        return x.muller.value_counts().idxmax()
    else:
        return "unknown_scaffold"


if __name__ == "__main__":
    DEBUG = False

    if DEBUG:
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("../../config/common.yaml")
        tag = config["tag"]
        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input={
                "orthologs": "../output/expression-atlas-wf/ortholog_annotation.feather",
                "gene_annotation": f"../references/gene_annotation_dmel_{tag}.feather",
                "primary2secondary": f"../references/primary2secondary_dmel_{tag}.pkl",
            },
        )
    main()
