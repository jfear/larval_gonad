import os

import re
import pandas as pd


def main():
    with open(snakemake.input[0]) as fin, open(snakemake.output[0], "w") as fout:
        for row in fin.read().strip().split("\n"):
            cell = snakemake.wildcards.sample + "_" + re.sub("-1", "", row) + "\n"
            fout.write(cell)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="adult-testis-scRNAseq-wf",
            input="../../output/adult-testis-scRNAseq-wf/Ral517/outs/raw_feature_bc_matrix/barcodes.tsv",
            wildcards=dict(sample="Ral517"),
        )

    main()
