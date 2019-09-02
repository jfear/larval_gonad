"""Creates a YOgn to Muller arm mapping.

Uses majority rule of D. melanogaster orthologs to assign scaffolds to
muller elements.

"""
import os
import pandas as pd
from larval_gonad.io import pickle_load, pickle_dump

MULLER = {"X": "A", "2L": "B", "2R": "C", "3L": "D", "3R": "E", "4": "F", "Y": "Y"}
UNK = "unknown_scaffold"


def main():
    # Muller arm assignment for dmel
    fbgn2muller = (
        pd.read_feather(snakemake.input.gene_annot)
        .set_index("FBgn")
        .FB_chrom.map(MULLER)
        .dropna()
        .rename("muller")
    )

    # Ortholog map from YO to updated dmel FBgn
    fbgn_updater = pickle_load(snakemake.input.primary2secondary)
    orthologs = {
        k: fbgn_updater[v]
        for k, v in pickle_load(snakemake.input.orthologs).items()
        if v in fbgn_updater
    }

    # Scaffold map to Muller element
    metadata = pd.read_feather(snakemake.input.yogn_annot)
    scaffold_to_muller = (
        metadata.loc[:, ["YOgn", "chrom"]]
        .assign(FBgn=lambda df: df.YOgn.map(orthologs))
        .query("FBgn == FBgn")
        .merge(fbgn2muller, on="FBgn")
        .groupby("chrom")
        .apply(consensus_muller)
        .to_dict()
    )

    # YOgn to Muller
    yogn2muller = metadata.set_index("YOgn").chrom.map(scaffold_to_muller).fillna(UNK).to_dict()

    pickle_dump(yogn2muller, snakemake.output[0])


def consensus_muller(x):
    if x.shape[0] > 20:
        return x.muller.value_counts().idxmax()
    else:
        return UNK


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input=dict(
                orthologs="../output/expression-atlas-wf/dana_YOgn_to_dmel_ortholog.pkl",
                yogn_annot="../output/expression-atlas-wf/dana_YOgn_metadata.feather",
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                primary2secondary="../references/primary2secondary_dmel_r6-26.pkl",
            ),
            wildcards=dict(species="dana"),
        )

    main()
