"""Playing with plots for galletta"""
#%%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from statsmodels.stats.weightstats import ttest_ind

plot_defaults = {
    'ylim': (-0.1, 4.3),
    'xlabel': 'Chromosome',
    'ylabel': 'Normalized Pixel Intensity\n(arbitrary unit)',
}

#%%
s2 = (
    pd.read_excel(
        "data/external/galletta/phos_over_tot_data_late_SC_only_061819.xlsx",
        sheet_name="S2",
        skiprows=1,
    )
    .melt(value_name="norm", var_name="chrom")
    .dropna()
    .assign(antibody="S2")
)

#%%
ax = sns.boxplot("chrom", "norm", data=s2)
ax.set(title='S2', **plot_defaults)
ttest_ind(s2.query('chrom == "X"').norm, s2.query('chrom == "A"').norm, usevar="unequal")

#%%
s5 = (
    pd.read_excel(
        "data/external/galletta/phos_over_tot_data_late_SC_only_061819.xlsx",
        sheet_name="S5",
        skiprows=1,
    )
    .melt(value_name="norm", var_name="chrom")
    .dropna()
    .assign(antibody="S5")
)

#%%
ax = sns.boxplot("chrom", "norm", data=s5)
ax.set(title='S5', **plot_defaults)
ttest_ind(s5.query('chrom == "X"').norm, s5.query('chrom == "A"').norm, usevar="unequal")
