import pandas as pd

from larval_gonad.io import feather_to_cluster_rep_matrix, melt_cluster_rep_matrix
from larval_gonad.normalization import zscore


TPM = snakemake.input[0]
OUTPUT_FILE = snakemake.output[0]


def main():
    tpm = feather_to_cluster_rep_matrix(TPM).dropna(axis=1, how='all')

    _zscore = zscore(tpm).dropna()

    melt_cluster_rep_matrix(_zscore, name='zscore').to_feather(OUTPUT_FILE)


if __name__ == "__main__":
    main()
