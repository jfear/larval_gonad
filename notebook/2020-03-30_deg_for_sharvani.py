# %%
import pickle
from pathlib import Path

import numpy as np
import pandas as pd


# %%
yo2fb = pickle.load(open("../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dmel.pkl", "rb"))
fbgn2symbol = pd.read_feather("../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()

# %%
TISSUES = [
    "AC",
    "DG",
    "GE",
    "GO",
    "HD",
    "RE",
    "TX",
    "WB"
]
# %%
for tissue in TISSUES:
    df = (
        pd.read_table(f"../output/expression-atlas-wf/sex_biased_expression/w1118_{tissue}.tsv")
        .assign(FBgn=lambda x: x.YOgn.map(lambda y: yo2fb.get(y, np.nan))).dropna().drop("YOgn", axis=1)
        .merge(fbgn2symbol, on="FBgn").set_index(["FBgn", "gene_symbol"])
    )
    df.to_csv(f"../output/expression-atlas-wf/sex_biased_expression_fb/w1118_{tissue}.tsv", sep="\t")

# %%
