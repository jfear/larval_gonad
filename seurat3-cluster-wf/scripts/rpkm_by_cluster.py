import pandas as pd

from larval_gonad.io import feather_to_cluster_matrix, melt_cluster_matrix
from larval_gonad.normalization import rpkm


RAW_AGG = snakemake.input["raw_agg"]
GENE_ANNOTATION = snakemake.input["gene_metadata"]

OUTPUT_FILE = snakemake.output[0]


def main():
    fbgn2len = (
        pd.read_feather(GENE_ANNOTATION, columns=["FBgn", "length"])
        .set_index("FBgn")
        .length
    )

    raw_agg = feather_to_cluster_matrix(RAW_AGG).dropna(axis=1, how="all")

    _rpkm = rpkm(raw_agg, fbgn2len).dropna()

    melt_cluster_matrix(_rpkm, name='RPKM').to_feather(OUTPUT_FILE)


if __name__ == "__main__":
    main()
