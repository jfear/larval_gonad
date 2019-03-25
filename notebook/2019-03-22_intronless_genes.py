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
# # Intronless genes

# %% [markdown]
# I want to make a list of genes that do not have introns.

# %%
import os
from pathlib import Path

from gffutils import FeatureDB

# %%
fname = Path(os.environ['REFERENCES_DIR'])  / 'dmel/r6-16/gtf/dmel_r6-16.gtf.db'
db = FeatureDB(fname.as_posix())

# %%
intron_counts = {}
for gene in db.features_of_type('gene'):
    intron_counts[gene.id] = 0

# %%
for intron in db.create_introns():
    intron_counts[intron.attributes['gene_id'][0]] += 1

# %%
for fbgn, cnt in intron_counts.items():
    if cnt == 0:
        print(fbgn)
        break

# %%
coding = set()
noncoding = set()
for gene in db.features_of_type('gene'):
    if len(list(db.children(gene, featuretype='CDS'))) == 0:
        noncoding.add(gene.id)
    else:
        coding.add(gene.id)

# %%
len(noncoding), len(coding)

# %%
with open('/home/fearjm/Downloads/FlyBase_IDs (1).txt', 'r') as fh:
    fbnc = fh.read().strip().split('\n')

# %%
fbnc = set(fbnc)

# %%

# %%
len(noncoding.intersection(fbnc))

# %%
fbnc - noncoding

# %%
