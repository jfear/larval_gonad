import os
import re

import numpy as np
import pandas as pd

from larval_gonad.io import cellranger_counts


def main():
    df = pd.concat([get_cell_metadata(fname) for fname in snakemake.input])
    df.rename_axis("cell_id").reset_index().to_feather(snakemake.output[0])


def get_cell_metadata(fname):
    rep_num = re.findall(r"testis(\d)", fname)[0]
    dat = cellranger_counts(fname)
    return pd.DataFrame(
        dict(
            nUMI=np.squeeze(np.asarray(dat.matrix.sum(axis=0))),
            nFeature=np.squeeze(np.asarray(np.sum(dat.matrix > 0, axis=0))),
        ),
        index=(f"rep{rep_num}_{x}" for x in dat.barcodes),
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="cellselection-wf",
            input=[
                "../output/cellranger3-wf/testis1/outs/raw_feature_bc_matrix.h5",
                "../output/cellranger3-wf/testis2/outs/raw_feature_bc_matrix.h5",
                "../output/cellranger3-wf/testis3/outs/raw_feature_bc_matrix.h5",
            ],
        )

    main()
