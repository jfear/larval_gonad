"""Calculate chromosome level stats for sex-biased expression."""
from typing import Tuple

import pandas as pd
from scipy.stats import chi2_contingency

from larval_gonad.stats import run_chisq

CHROM = snakemake.params.chrom_order


def main():
    df = (
        pd.read_feather(snakemake.input[0])
        .groupby(["species", "tissue"])
        .apply(cross_tabulate)
    )

    full_stats = (
        df.groupby(["species", "tissue"], as_index=False).apply(run_chisq).droplevel(0)
    )
    full_stats.to_csv(snakemake.output.full_stats, sep="\t")

    filtered_stats = df.groupby(["species", "tissue"], as_index=False).apply(
        filter_stats
    )
    filtered_stats.to_csv(snakemake.output.filtered_stats, sep="\t")


def cross_tabulate(df: pd.DataFrame) -> pd.DataFrame:
    # Build Crosstab table
    chrom = df.chrom
    bias = df.bias
    ct = pd.crosstab(chrom, bias).fillna(0)
    ct.columns.name = ""
    return ct


def filter_stats(ct: pd.DataFrame) -> Tuple[float, pd.DataFrame]:
    """Calculates chi-square and post hoc tests."""
    # Overall table chi-square
    _, table_pval, _, _ = chi2_contingency(ct)

    # Post Hoc Tests
    res = (
        run_chisq(ct)
        .loc[(slice(None), slice(None), slice(None), "fdr q-value"), :]
        .droplevel(-1)
    )
    res.columns = [f"qval_{x}" for x in res.columns]

    # Return Table and Post Hoc tests
    res["table_pval"] = table_pval
    return res


if __name__ == "__main__":
    main()
