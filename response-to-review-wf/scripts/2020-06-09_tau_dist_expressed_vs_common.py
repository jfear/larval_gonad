#%%
import joblib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%
expressed = joblib.load("../output/cellselection-wf/expressed_genes.pkl")
common = joblib.load("../output/cellselection-wf/commonly_expressed_genes.pkl")

# %%
tau = pd.read_feather("../output//expression-atlas-wf/dmel_tau.feather").set_index("FBgn").male_tau
tsps = pd.read_feather("../output//expression-atlas-wf/dmel_tsps.feather").set_index("FBgn").male_tsps

# %%
sns.kdeplot(tau[expressed], label="All Expressed")
sns.kdeplot(tau[common], label="Commonly Expressed")
plt.suptitle("tau Score")



# %%
sns.kdeplot(tsps[expressed], label="All Expressed")
sns.kdeplot(tsps[common], label="Commonly Expressed")
plt.suptitle("TSPS Score")
