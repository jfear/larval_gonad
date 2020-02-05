import os
import pandas as pd

from larval_gonad.normalization import tpm


def main():
    df = pd.concat(
        [read_l3_sc(), read_l3_bulk()], sort=True, axis=1, join="inner"
    )  # type: pd.DataFrame

    df.corr(method="spearman").rename_axis("Spearman rho").to_csv(snakemake.output[0], sep="\t")


def read_l3_sc() -> pd.DataFrame:
    gene_lengths = pd.read_feather(snakemake.input.gene_annot).set_index("FBgn")["length"].squeeze()

    raw = (
        pd.read_feather(snakemake.input.larval_scrnaseq)
        .melt(id_vars="FBgn", var_name="cell_id", value_name="UMI")
        .assign(rep=lambda x: x.cell_id.str.extract(r"(rep\d)_.*", expand=False))
        .groupby(["FBgn", "rep"]).UMI.sum()
        .query("rep != 'rep4")
    )

    raw_wide = pd.pivot_tale(raw, values="UMI", index="FBgn", columns="rep")

    raw_norm = tpm(raw_wide, gene_lengths).dropna()
    raw_norm.columns = [f"l3_scRNAseq_{x}" for x in raw_norm_columns]

    return raw_norm


def read_l3_bulk() -> pd.DataFrame:
    df = pd.read_csv(snakemake.input.larval_bulk, sep="\t", index_col=0).rename_axis("FBgn")
    cols = [f"l3_bulk_{x}" for x in df.columns]
    df.columns = cols
    return df


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
                larval_scrnaseq="../../output/cellselection-wf/raw.feather",
                larval_bulk="../../output/bulk2-rnaseq-wf/rnaseq_aggregation/tpm_gene_level_counts.tsv",
            ),
        )

    main()
