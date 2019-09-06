import os
import numpy as np
import pandas as pd

from svgutils.compose import Figure, Panel, SVG

from larval_gonad.io import safe_gene_name


def main():
    fbgn2symbol = (
        pd.read_feather(snakemake.input[0], columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .squeeze()
    )

    # Get gene symbols
    pairs = [(safe_gene_name(fbgn2symbol[fbgn]), fbgn) for fbgn in snakemake.params.genes]

    # Calculate the nubmer of panels
    ncol = 4
    nrow = np.ceil(len(pairs) / 4)

    width = "21cm"
    height = f"{nrow * 5}cm"

    Figure(
        width,
        height,
        *(SVG(snakemake.params.pattern.format(symbol=x[0], FBgn=x[1])).scale(0.5) for x in pairs),
    ).tile(ncol, nrow).save(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        tag = config["tag"]
        lit_genes = read_config("config/literature_genes.yaml")

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=f"../references/gene_annotation_dmel_{tag}.feather",
            params=dict(
                pattern="../output/seurat3-cluster-wf/combined_n3_figures/gene_projections/{symbol}_{FBgn}.svg",
                genes=lit_genes["PS"],
            ),
            wildcards=dict(group="PS"),
        )

    main()
