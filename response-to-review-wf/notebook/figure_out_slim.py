# %% [markdown]
# # Figure out custom SLIM
#
# BO wants to make a summary ribbon like on FlyBase. I order to do that I need
# a custom slim based on the Drosophila obo.

# %%
from goatools.godag_plot import plot_gos
from goatools.mapslim import mapslim
from larval_gonad.gene_ontology import slimDag, oboDag

# %% [markdown]
# ## Molecular Function

# %%
fb_mf = [
    "GO:0003824",
    "GO:0098772",
    "GO:0038023",
    "GO:0005102",
    "GO:0005215",
    "GO:0005198",
    "GO:0008092",
    "GO:0003677",
    "GO:0003723",
    "GO:0140110",
    "GO:0036094",
    "GO:0046872",
    "GO:0008289",
    "GO:0030246",
]

fb_bp = [
    "GO:0008283",
    "GO:0007049",
    "GO:0071840",
    "GO:0051234",
    "GO:0033036",
    "GO:0032502",
    "GO:0000003",
    "GO:0002376",
    "GO:0050877",
    "GO:0007610",
    "GO:0050896",
    "GO:0023052",
    "GO:0010467",
    "GO:0019538",
    "GO:0006259",
    "GO:0044281",
]

fb_cc = [
    "GO:0005576",
    "GO:0005829",
    "GO:0005856",
    "GO:0005739",
    "GO:0005634",
    "GO:0005694",
    "GO:0016020",
    "GO:0071944",
    "GO:0042995",
    "GO:0005773",
    "GO:0012505",
    "GO:0045202",
    "GO:0032991",
]

fb = set(fb_bp + fb_cc + fb_mf)

# %%
def get_parents(term) -> set:
    parents = set()
    paths = oboDag.paths_to_top(term)
    for path in paths:
        for term_ in path:
            if term_.id == term:
                continue
            parents.add(term_.id)
    return parents


def run(ribbons: set):
    slimSet = set(slimDag)
    to_keep = ribbons.intersection(slimSet)
    to_add = ribbons - slimSet

    to_remove = set()
    for slim in slimSet - to_keep:
        parents = get_parents(slim)
        if parents.intersection(to_add):
            to_remove.add(slim)

    for slim in slimSet - to_remove - to_add:
        to_keep.add(slim)

    return to_keep, to_add, to_remove

to_keep, to_add, to_remove = run(fb)

# %%
print(len(to_keep))
print(len(to_add))
print(len(to_remove))

# %%
to_keep.intersection(to_add), to_keep.intersection(to_remove), to_remove.intersection(to_add)

# %%
sorted(to_remove)
