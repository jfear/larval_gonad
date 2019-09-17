import pandas as pd


def get_fbgn2chrom(gene_metadata: str) -> pd.Series:
    return (
        pd.read_feather(gene_metadata, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .rename("chrom")
    )


def get_fbgn2symbol(gene_metadata: str) -> pd.Series:
    return (
        pd.read_feather(gene_metadata, columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()
    )


def get_fbgn2length(gene_metadata: str) -> pd.Series:
    return pd.read_feather(gene_metadata, columns=["FBgn", "length"]).set_index("FBgn").squeeze()

