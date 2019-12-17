import os
from pathlib import Path

import pandas as pd
import scipy.io


def main():
    barcodes = pd.read_csv(snakemake.input.barcodes, sep="\t", header=None, index_col=0).index
    barcodes.name = "cell_id"

    features = pd.read_csv(snakemake.input.features, sep="\t", header=None, index_col=0).index
    features.name = "FBgn"

    matrix = scipy.io.mmread(snakemake.input.matrix)

    df = pd.DataFrame(matrix.todense(), index=features, columns=barcodes)

    df.reset_index().to_feather(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="adult-testis-scRNAseq-wf",
            input=dict(
                barcodes="../output/adult-testis-scRNAseq-wf/Ral517/barcodes.tsv",
                features="../output/adult-testis-scRNAseq-wf/Ral517/features.tsv",
                matrix="../output/adult-testis-scRNAseq-wf/Ral517/matrix.mtx",
            ),
        )

    main()
