"""Create a temporary mapping of dmel FBgn to Muller.

This script converts the YOgn to Muller element and D. melanogaster FBGn to
Muller element. It uses the YOgn to D. mel FBgn mapping to update the YOgn
to Muller element mapping.

"""
import os
import pandas as pd

from larval_gonad.io import pickle_load, pickle_dump

def main():
    muller = pickle_load(snakemake.input.muller)
    ortholog = pickle_load(snakemake.input.ortholog)

    dat = pd.Series({
        ortholog[k]: v
        for k, v in muller.items()
        if k in ortholog
    }, name=snakemake.wildcards.species)
    dat.index.name = "FBgn"
    dat.to_pickle(snakemake.output[0])


if __name__ == '__main__':
    if os.getenv('SNAKE_DEBUG', False):
        from larval_gonad.debug import snakemake_debug
        snakemake = snakemake_debug(
            workdir='expression-atlas-wf',
            input=dict(
                muller="../output/expression-atlas-wf/YOgn_to_muller/dmel.pkl",
                ortholog="../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dmel.pkl"
            ),
            wildcards=dict(species="dmel")
        )

    main()
