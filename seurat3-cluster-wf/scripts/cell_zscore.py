import os
import pandas as pd

from larval_gonad.normalization import zscore


def main():
    norm = pd.read_feather(snakemake.input[0]).set_index("FBgn").drop("gene_symbol", axis=1)
    zscore(norm).reset_index().to_feather(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input="../output/seurat3-cluster-wf/combined_n3_normalized.feather",
        )

    main()
