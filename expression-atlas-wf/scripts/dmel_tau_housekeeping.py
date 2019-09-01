import os
from functools import partial
import pandas as pd

from larval_gonad.io import pickle_load, pickle_dump


def main():
    # Load mapping of YOgn to FBgn
    annot = pickle_load(snakemake.input.annot[0])

    pickle_dump(intersect_fbgns(snakemake.input.male, annot), snakemake.output.male)
    pickle_dump(intersect_fbgns(snakemake.input.female, annot), snakemake.output.female)


def intersect_fbgns(file_names, annot):
    return list(set.intersection(*list(map(partial(convert_to_fbgn, annot=annot), file_names))))


def convert_to_fbgn(file_name, annot):
    return set(
        [
            fbgn
            for fbgn in map(lambda x: annot.get(x, None), pickle_load(file_name))
            if fbgn is not None
        ]
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input=dict(
                male=[
                    "../output/expression-atlas-wf/tau_housekeeping/w1118_male.pkl",
                    "../output/expression-atlas-wf/tau_housekeeping/orgR_male.pkl",
                ],
                female=[
                    "../output/expression-atlas-wf/tau_housekeeping/w1118_female.pkl",
                    "../output/expression-atlas-wf/tau_housekeeping/orgR_female.pkl",
                ],
                annot="../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dmel.pkl",
            ),
        )

    main()
