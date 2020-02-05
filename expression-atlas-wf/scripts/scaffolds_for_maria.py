import os
from pathlib import Path
import pickle

import numpy as np
import pandas as pd


def main():
    orthologs = {Path(x).stem: x for x in snakemake.input.orthologs}
    scaffolds = {Path(x).stem: x for x in snakemake.input.scaffolds}

    # Merge dmel seprately to avoid including non-orhtologs.
    dmel = read_data("dmel", orthologs["dmel"], scaffolds["dmel"])

    other_species = pd.concat(
        [
            read_data(species, ortholog, scaffolds[species])
            for species, ortholog in orthologs.items()
            if species != "dmel"
        ],
        axis=1,
        sort=False,
    ).rename_axis("dmel_FBgn")

    df = other_species.join(dmel, how="inner")
    df.to_csv(snakemake.output[0], sep="\t")


def read_data(species, ortholog, scaffold):
    mapper = pickle.load(open(ortholog, "rb"))

    def map_fbgn(x):
        return mapper.get(x, np.nan)

    df = (
        pd.read_feather(scaffold)
        .assign(dmel_FBgn=lambda x: x.YOgn.map(map_fbgn))
        .dropna()
        .set_index("dmel_FBgn")
        .drop("YOgn", axis=1)
    )
    df.columns = [f"{species}_{x}" for x in df.columns]
    return df


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                orthologs=[
                    "../../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dana.pkl",
                    "../../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dmel.pkl",
                ],
                scaffolds=[
                    "../../output/expression-atlas-wf/YOgn_metadata/dana.feather",
                    "../../output/expression-atlas-wf/YOgn_metadata/dmel.feather",
                ],
            )
        )

    main()
