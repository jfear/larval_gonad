"""Create mapping dictionary from YOgn to dmel FBgn.

Reads the ortholog table and creates a dictionary mapping from YOgn to dmel's
FBgn.

"""
import os
import pandas as pd

from larval_gonad.io import pickle_dump


def main():
    yo_to_ortho = (
        pd.read_csv(snakemake.input[0], sep="\t")
        .pipe(lambda x: x[x.YOgnID != "-"])
        .set_index("YOgnID")
        .Dmel.rename("FBgn")
        .to_dict()
    )

    pickle_dump(yo_to_ortho, snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input="../output/expression-atlas-wf/orthologs/dmel.tsv",
            wildcards=dict(species="dmel"),
        )

    main()
