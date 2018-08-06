"""Figure 1. Identification, annotation and validation of clusters."""

from svgutils.compose import *
from common import label

Figure(
    #"21cm", "10.16cm",
    #"792", "612",
    "576", "288",
#     Panel(SVG(snakemake.input.dia).scale(1), label("a")),
    Panel(SVG(snakemake.input.tsne).scale(1), label("b")).move(0, 144),
    Panel(SVG(snakemake.input.hmap_all).scale(1), label("c")).move(144, 0),
    Panel(SVG(snakemake.input.hmap_lit).scale(1), label("d")).move(288, 0),
    Panel(SVG(snakemake.input.hmap_ptrap).scale(1), label("e")).move(288, 144),
#     Grid(10, 10)
).save(snakemake.output[0])
