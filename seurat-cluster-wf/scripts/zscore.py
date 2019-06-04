"""Calculate Z-scores

Calculates the row based z-score. This is useful for making heatmaps.
"""
from pathlib import Path

import pandas as pd

from larval_gonad.normalization import zscore


def main():
    df = pd.read_feather(snakemake.input[0]).set_index('FBgn')
    zscore(df).to_csv(snakemake.output[0], sep='\t')


if __name__ == "__main__":
    main()
