import pandas as pd
import joblib

ARMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]


def _debug():
    from larval_gonad.mock import MockSnake

    snakemake = MockSnake(
        input=dict(
            annotation="../../references/gene_annotation_dmel_r6-26.feather",
            expressed="../../output/cellselection-wf/expressed_genes.pkl",
            widely_expressed="../../output/cellselection-wf/commonly_expressed_genes.pkl",
        )
    )


def main():
    annotation = pd.read_feather(snakemake.input.annotation).query("FB_chrom in @ARMS")
    expressed = joblib.load(snakemake.input.expressed)
    widely_expressed = joblib.load(snakemake.input.widely_expressed)

    gene_counts = pd.concat([
        annotation.groupby("FB_chrom").size().rename("Number of Genes (FB6.26)"),
        annotation.query("FBgn in @expressed").groupby("FB_chrom").size().rename("Number of Expressed Genes"),
        annotation.query("FBgn in @widely_expressed").groupby("FB_chrom").size().rename("Number of Widely Expressed Genes"),
    ], axis=1).fillna(0).astype(int).rename_axis("Chromosome Arm")

    gene_counts.to_csv(snakemake.output[0], sep="\t")





if __name__ == "__main__":
    main()
