import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

pd.set_option("display.float_format", "{:.6f}".format)


def main():

    qvals = (
        pd.read_feather(snakemake.input[0])
        .groupby(["stage", "strain", "tissue", "chrom"])
        .apply(run_ttest, "less")
        .dropna()
        .pipe(adjust_pvalues)
        .to_frame()
        .unstack()
        .fillna(1)
        .astype(float)
        .droplevel(level=0, axis=1)
        .loc[:, ["X", "2L", "2R", "3L", "3R", "4", "Y"]]
    )

    qvals.to_csv(snakemake.output[0], sep="\t")


def run_ttest(df: pd.DataFrame, alternative="less") -> float:
    male = df.query("sex == 'm'")["Avg TPM Per Chromosome"].values
    female = df.query("sex == 'f'")["Avg TPM Per Chromosome"].values
    stat, pval = ttest_ind(male, female, equal_var=False)

    if alternative == "less":
        if stat < 0:
            return pval / 2
        else:
            return 1 - (pval / 2)
    elif alternative == "greater":
        if stat > 0:
            return pval / 2
        else:
            return 1 - (pval / 2)
    return pval


def adjust_pvalues(sr: pd.Series) -> pd.Series:
    _, qvals, _, _ = multipletests(sr, method="fdr_bh")
    return pd.Series(qvals, index=sr.index, name="q-value")


if __name__ == "__main__":
    main()
