from pathlib import Path

import pandas as pd
import scipy.io


INPUT_FILES = snakemake.input
OUTPUT_FILE = snakemake.output[0]


def main():
    df = pd.concat((create_matrix(mtx) for mtx in INPUT_FILES), axis=1)
    df.reset_index().to_feather(OUTPUT_FILE)


def create_matrix(mtx):
    genes = pd.read_csv(
        Path(mtx).parent / "genes.tsv", sep="\t", header=None, index_col=0
    ).index
    genes.name = "FBgn"

    barcodes = pd.read_csv(
        Path(mtx).parent / "barcodes.tsv", sep="\t", header=None, index_col=0
    ).index
    barcodes.name = "cell_id"

    matrix = scipy.io.mmread(mtx)

    return pd.DataFrame(matrix.todense(), index=genes, columns=barcodes)


if __name__ == "__main__":
    main()
