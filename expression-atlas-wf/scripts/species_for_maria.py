import os
from pathlib import Path

import pandas as pd

from larval_gonad.io import pickle_load


def main():
    df = pd.concat([muller_table(), fbgn_to_yogn()], axis=1, sort=True)
    df.to_csv(snakemake.output[0], sep="\t")


def muller_table():
    return (
        pd.read_feather(snakemake.input.muller)
        .set_index("FBgn")
        .rename_axis("dmel_FBgn")
        .rename(columns=lambda x: f"{x}_muller")
    )


def fbgn_to_yogn():
    return pd.concat(read_pickels(), axis=1, sort=True).rename_axis("dmel_FBgn")


def read_pickels():
    for species in snakemake.input.yogns:
        name = Path(species).stem
        dat = pickle_load(species)
        yield pd.Series({v: k for k, v in dat.items()}, name=f"{name}_YOgn")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                muller="../../output/expression-atlas-wf/aggregated_muller_table.feather",
                yogns=[
                    "../../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dana.pkl",
                    "../../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dyak.pkl",
                ],
            )
        )

    main()
