"""Figure 1. Identification, annotation and validation of clusters."""

from svgutils.compose import *
from common import label

Figure(
    #"21cm", "10.16cm",
    #"792", "612",
    "432", "288",
#     Panel(SVG(snakemake.input.dia), label("A")),
    Panel(
        SVG(snakemake.input.tsne),
        label("B")
    ).move(0, 144),
    Panel(
        SVG(snakemake.input.hmap_all),
        label("C")
    ).move(144, 0),
    Panel(
        SVG(snakemake.input.hmap_lit),
        label("D")
    ).move(288, 0),
    Panel(
        SVG(snakemake.input.hmap_ptrap),
        label("E")
    ).move(288, 144),
#     Grid(10, 10)
).save(snakemake.output[0])
