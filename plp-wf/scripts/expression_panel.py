"""Trial gene panel.

Brian G. has a list of genes that he is focused on, I wanted to pull them out
and make a nice panel for his supplement.

"""
import pandas as pd
from svgutils.compose import Panel, SVG, Figure
from larval_gonad.io import safe_gene_name

fbgn2symbol = (
    pd.read_feather("references/gene_annotation_dmel_r6-26.feather").set_index("FBgn").gene_symbol
)

trial_fbgns = pd.read_csv("data/external/galletta/trial_list.txt", sep="\t").FBGN.tolist()

trial_genes = fbgn2symbol.reindex(trial_fbgns)


Figure(
    "12cm",
    "12cm",
    *[
        SVG(
            f"output/seurat3-cluster-wf/combined_n3_figures/expression_patterns/{safe_gene_name(symbol)}_{fbgn}.svg"
        ).scale(0.2)
        for fbgn, symbol in trial_genes.iteritems()
    ],
).tile(5, 11).save("output/plp-wf/expression_panel.svg")


Figure(
    "12cm",
    "12cm",
    *[
        SVG(
            f"output/seurat3-cluster-wf/combined_n3_figures/gene_projections/{safe_gene_name(symbol)}_{fbgn}.svg"
        ).scale(0.1)
        for fbgn, symbol in trial_genes.iteritems()
    ],
).tile(5, 11).save("output/plp-wf/projection_panel.svg")
