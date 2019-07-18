"""Determine chromosome consensus for scaffolds.

We want to include the out groups D. willistoni, D. virilis, and D.
mojavensus in our analysis. However, genes are assigned to scaffolds and not
Muller elements. We assign chromosome arm based on the consensus from D.
melanogaster.

"""
import numpy as np
import pandas as pd

ORTHOLOGS = snakemake.input["orthologs"]
GENE_ANNOTATION = snakemake.input["gene_annotation"]
OUTPUT_FILE = snakemake.output[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), "neox-wf/scripts"))
#     print(os.getcwd())
# except:
#     pass
# ORTHOLOGS = '../../output/neox-wf/ortholog_annotation.feather'
# GENE_ANNOTATION = '../../references/gene_annotation_dmel_r6-24.feather'

DMEL_MULLER = {"X": "A", "2L": "B", "2R": "C", "3L": "D", "3R": "E", "4": "F"}
DPSE_MULLER = {"XL": "A", "4": "B", "3": "C", "XR": "D", "2": "E", "5": "F"}


def main():
    fbgn2muller = (
        pd.read_feather(GENE_ANNOTATION, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .map(DMEL_MULLER)
        .rename("muller")
    )

    orthologs = pd.read_feather(ORTHOLOGS).set_index("FBgn").join(fbgn2muller)

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
