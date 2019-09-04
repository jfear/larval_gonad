import os
import pandas as pd

import matplotlib.pyplot as plt
from svgutils.compose import Figure, Panel, SVG, Text


def main():
    label_defaults = dict(x=35, y=10, size=12, weight="bold", font="Helvetica")
    annot_defaults = label_defaults.copy()
    annot_defaults.update(dict(x=520, y=80, size=10))

    Figure(
        "20.32cm",
        "13cm",
        Panel(
            SVG(snakemake.input.wb),
            Text("A", **label_defaults),
            Text("Whole Body", **annot_defaults),
        ),
        Panel(
            SVG(snakemake.input.go), Text("B", **label_defaults), Text("Gonad", **annot_defaults)
        ),
        Panel(
            SVG(snakemake.input.ac), Text("C", **label_defaults), Text("Carcass", **annot_defaults)
        ),
    ).tile(1, 3).save(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                wb="../output/science_submission/demas/WB_x_and_4th.svg",
                go="../output/science_submission/demas/GO_x_and_4th.svg",
                ac="../output/science_submission/demas/AC_x_and_4th.svg",
            ),
            output="../output/science_submission/demas/test.svg",
        )

    main()
