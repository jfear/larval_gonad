"""Figure 1. Identification, annotation and validation of clusters."""

from svgutils.compose import *
from common import label

Figure(
    "21cm", "10.16cm",
    Panel(SVG(snakemake.input.tsne).scale(1.6), label("b")).move(0, 175),
    Panel(SVG(snakemake.input.hmap_all).scale(1.7), label("c")).move(200, 0),
    Panel(SVG(snakemake.input.hmap_lit).scale(1.6), label("d")).move(450, 0),
    #Panel(SVG(snakemake.input.comp).scale(0.9).move(10, 0), label("c")).move(560, 0),
    #Panel(SVG(snakemake.input.age).scale(0.9).move(19, 0)).move(560, 90),
     Grid(20, 20)
).save(snakemake.output[0])
