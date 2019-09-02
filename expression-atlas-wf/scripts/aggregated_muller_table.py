"""Creates single table with all species Muller elements"""
import os
import pandas as pd

from larval_gonad.io import pickle_load


def main():
    df = pd.concat((pickle_load(file_name) for file_name in snakemake.input), axis=1, sort=True)
    df.index.name = "FBgn"
    df.reindex(snakemake.params.species, axis="columns").reset_index().to_feather(
        snakemake.output[0]
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input=[
                "../output/expression-atlas-wf/FBgn_to_muller_dmel.pkl",
                "../output/expression-atlas-wf/FBgn_to_muller_dana.pkl",
            ],
            params=dict(species=["dmel", "dana"]),
        )

    main()
