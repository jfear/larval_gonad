"""Determine chromosome consensus for scaffolds.

We want to include the out groups D. willistoni, D. virilis, and D.
mojavensus in our analysis. However, genes are assigned to scaffolds and not
Muller elements. We assign chromosome arm based on the consensus from D.
melanogaster.

"""
import pickle
import numpy as np
import pandas as pd

ORTHOLOGS = snakemake.input["orthologs"]
GENE_ANNOTATION = snakemake.input["gene_annotation"]
PRIMARY_SECONDARY = snakemake.input['primary2secondary']
OUTPUT_FILE = snakemake.output[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), "expression-atlas-wf/scripts"))
#     print(os.getcwd())
# except:
#     pass
# ORTHOLOGS = '../../output/expression-atlas-wf/ortholog_annotation.feather'
# GENE_ANNOTATION = '../../references/gene_annotation_dmel_r6-26.feather'
# PRIMARY_SECONDARY = '../../references/primary2secondary_dmel_r6-26.pkl'

DMEL_MULLER = {"X": "A", "2L": "B", "2R": "C", "3L": "D", "3R": "E", "4": "F"}
DPSE_MULLER = {"XL": "A", "4": "B", "3": "C", "XR": "D", "2": "E", "5": "F"}


def main():
    fbgn2muller = (
        pd.read_feather(GENE_ANNOTATION, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .map(lambda x: DMEL_MULLER.get(x, 'unknown_scaffold'))
        .rename("muller")
    )

    with open(PRIMARY_SECONDARY, 'rb') as fh:
        primary2secondary = pickle.load(fh)

    orthologs = pd.read_feather(ORTHOLOGS).set_index("FBgn")

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

    dpse_sanity = pd.concat([dpse_fbgn2annotated, dpse_fbgn2consensus], axis=1).dropna(subset=['annotated'])

    # Assure that dpse annotation and consensus are identical
    assert (dpse_sanity.annotated.fillna("NAN") == dpse_sanity.consensus.fillna("NAN")).all()

    # Get muller for other samples
    dvir_muller = orthologs.dvir.map(
        orthologs.groupby("dvir").apply(consensus_muller).rename("muller")
    )
    dmoj_muller = orthologs.dmoj.map(
        orthologs.groupby("dmoj").apply(consensus_muller).rename("muller")
    )
    dwil_muller = orthologs.dwil.map(
        orthologs.groupby("dwil").apply(consensus_muller).rename("muller")
    )

    df_muller = pd.concat(
        [
            orthologs.muller.rename("dmel"),
            dpse_fbgn2consensus.rename("dpse"),
            dvir_muller,
            dmoj_muller,
            dwil_muller,
        ],
        axis=1,
    )

    df_muller.reset_index().to_feather(OUTPUT_FILE)


def consensus_muller(x):
    if x.shape[0] > 20:
        return x.muller.value_counts().idxmax()
    else:
        return 'unknown_scaffold'


if __name__ == "__main__":
    main()
