# %%
import pandas as pd

# %%
import os
try:
    os.chdir(os.path.join(os.getcwd(), 'notebook'))
    print(os.getcwd())
except:
    pass


# %%
df = pd.read_csv("../output/science_submission/cell_metadata.tsv", sep="\t", index_col=0)
df

# %%
df.is_cell.sum()
# %%
(df[df.is_cell ^ df.scrublet_is_multi].nFeature > 5000).sum()

# %%
(df[df.is_cell].scrublet_is_multi).sum()

# %%
df.flag_cell_used_in_study.sum()
