import os
import numpy as np
import pandas as pd

from svgutils.compose import Figure, Panel, SVG


def main():
    inputs = list(snakemake.input)

    # Calculate the nubmer of panels
    ncol = 4
    nrow = np.ceil(len(inputs) / 4)

    width = "21cm"
    height = f"{nrow * 5}cm"

    Figure(
        width,
        height,
        *(SVG(file_name).scale(0.5) for file_name in inputs),
    ).tile(ncol, nrow).save(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        tag = config["tag"]
        lit_genes = read_config("config/literature_genes.yaml")

        snakemake = snakemake_debug(
            workdir="germ-trajectory-wf",
            input=[
                "../output/germ-trajectory-wf/gene_projectsion/bol_FBgn0011206.svg",
                "../output/germ-trajectory-wf/gene_projectsion/tj_FBgn0000964.svg",
            ],
        )

    main()
