# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.0
#   kernelspec:
#     display_name: Python [conda env:larval_gonad]
#     language: python
#     name: conda-env-larval_gonad-py
# ---

# %% [markdown]
# Cleaning up protein trap scores from miriam.

# %%
from pickle import load
from io import StringIO
import re
import pandas as pd

# %%
with open('../output/science_submission/symbol2fbgn.pkl', 'rb') as fh:
    symbol2fbgn = load(fh)

# %%
miriam_scores = """gene_symbol,SP,ES,MS,LS,C,EC,MC,LC,PC,TE
ADD1,1,2,3,3,1,0,0,0,0,1
Ance,1,2,2,2,1,2,2,2,1,1
bol,1,1,2,3,0,0,0,0,0,0
CG17646,0,0,0,0,1,1,1,1,1,1
CG9747,0,1,1,1,0,1,1,1,1,1
Cht5,0,0,0,0,0,1,1,2,0,0
cindr,1,1,1,1,1,1,1,1,2,2
cmpy,0,0,0,0,0,1,1,1,0,0
Dek,2,2,2,2,1,1,1,1,2,2
dpr17,1,1,1,0,0,0,0,0,0,0
e(y)3,2,2,2,2,1,1,1,1,2,2
DIP-delta,2,2,1,1,0,0,0,0,1,2
Efa6,1,1,1,1,1,2,2,2,1,1
Efa6,2,1,1,1,2,2,2,2,1,1
Fas3,1,1,1,1,1,1,1,1,1,3
Fas3,0,0,0,0,0,0,0,0,0,2
fln,1,1,1,1,2,2,2,2,1,1
foxo,2,2,2,2,1,1,1,1,1,2
Fs(2)Ket,2,2,2,2,1,1,1,1,1,1
Mapmodulin,2,2,2,1,1,1,1,1,2,2
mbl,1,2,2,2,1,0,0,0,1,3
Nlg3,1,1,1,1,1,2,2,1,0,1
Nlg3,0,0,0,0,0,1,1,1,0,0
nord,1,1,1,1,0,0,0,0,0,2
Nrg,2,1,1,1,2,2,2,2,2,2
osa,1,1,1,1,2,2,2,2,2,2
p53,2,2,1,0,0,0,0,0,0,1
Pdcd4,3,3,3,3,1,2,2,2,3,3
Piezo,0,0,0,0,0,0,0,0,0,3
rdo,1,1,1,1,2,2,2,2,1,1
rdo,1,1,1,1,2,2,3,3,1,1
rdo,1,1,1,1,2,2,3,3,1,1
RunxA,1,1,1,1,1,1,1,1,1,1
Sap-r,1,1,1,1,2,3,3,3,2,3
sosie,1,1,1,1,1,2,2,2,1,1
spir,1,1,1,1,0,0,0,0,0,0
SRPK,2,2,2,2,0,0,0,0,1,1
stai,3,2,2,2,2,2,2,2,1,3
Syn,1,1,1,1,1,1,1,2,1,1
Syn,1,1,1,1,1,1,1,1,1,1
Tep2,0,1,1,1,0,2,2,2,2,0
tok,1,1,1,1,1,0,0,0,1,2
tutl,1,1,1,1,0,0,0,0,1,1
twin,1,1,1,1,0,0,0,0,0,0
"""

# %%
df = (
    pd.read_csv(StringIO(miriam_scores))
    .assign(FBgn=lambda df: df.gene_symbol.map(symbol2fbgn))
    .set_index('FBgn')
    .drop('gene_symbol', axis=1)
)

# %%
print(df.to_string())

# %%
