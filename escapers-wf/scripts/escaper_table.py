"""Pull out genes that go up in L1 compared to G."""
import pandas as pd


def main():
    deg = pd.read_feather(snakemake.input.deg).pipe(add_gene_annot).pipe(add_direction)
    deg.query("chrom == 'Y'").sort_values("FBgn")
    deg.query("chrom == '4'").sort_values("FBgn")

    pd.crosstab(deg.direction, deg.chrom)[snakemake.params.chrom_order]

    pd.options.display.max_rows = 100

    raw = (
        pd.read_feather("../../output/seurat3-cluster-wf/raw_by_cluster.feather")
        .pipe(add_gene_annot)
        .pipe(pivot)
        .pipe(get_germline)
        .pipe(filter_low_expressers)
        .reset_index()
        .sort_values("FBgn")
        .set_index(["chrom", "FBgn", "gene_symbol"])
    )
    raw

    raw.query("chrom == 'X'")
    raw.query("chrom == 'Y'")
    raw.query("chrom == '4'")
    raw.query("chrom == '4' and G > LPS")
    raw.query("chrom == '4' and G < LPS")


def add_gene_annot(df: pd.DataFrame) -> pd.DataFrame:
    gene_annot = pd.read_feather(
        snakemake.input.gene_annot, columns=["FBgn", "FB_chrom", "gene_symbol"]
    ).rename(columns={"FB_chrom": "chrom"})
    return df.merge(gene_annot)


def add_direction(df: pd.DataFrame, cutoff: float = 0.01) -> pd.DataFrame:
    up_in_gonia = (df.p_val_adj <= cutoff) & (df.avg_logFC > 0)
    up_in_late = (df.p_val_adj <= cutoff) & (df.avg_logFC < 0)
    df["direction"] = "NS"
    df.loc[up_in_gonia, "direction"] = "Gonia-Biased"
    df.loc[up_in_late, "direction"] = "Late Spermatocyte-Biased"
    return df


def pivot(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.set_index(["chrom", "FBgn", "gene_symbol", "cluster"])
        .unstack()
        .droplevel(0, axis=1)
    )


def get_germline(df: pd.DataFrame) -> pd.DataFrame:
    _df = df[["G", "EPS", "MPS", "LPS"]].copy()
    _df.columns = [col for col in _df.columns]
    return _df


def filter_low_expressers(df: pd.DataFrame, cutoff: int = 10):
    mask = df.sum(axis=1) >= cutoff
    print(mask.shape)
    return df[mask]


if __name__ == "__main__":
    main()
