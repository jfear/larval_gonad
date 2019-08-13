"""Creates a shelve for easier plotting.

Creates a data set ready for plotting.

"""
import shelve

import numpy as np
import pandas as pd


def main(snake: dict):
    ratios = munge_ratios(snake["ratios"])

    x_locations = get_cluster_codes(ratios.cluster)
    y_locations = (
        ratios.groupby(["cluster", "ratio_type"])
        .ratio.apply(get_whisker_location)
        .rename("y")
        .reset_index()
    )

    pvalues = munge_pvalues(snake["pvalues"], x_locations, y_locations)

    with shelve.open(snake["output_prefix"].replace(".dat", "")) as db:
        db["data"] = ratios
        db["pvalues"] = pvalues


def munge_ratios(input_file: str) -> pd.DataFrame:
    """Melt different ratios types into stacked format."""
    return (
        pd.read_feather(input_file)
        .melt(id_vars=["cell_id", "cluster", "rep"], var_name="ratio_type", value_name="ratio")
        .set_index("cell_id")
    )


def get_cluster_codes(cluster: pd.Categorical) -> pd.Series:
    """Get the X location for plotting p-value string."""
    categories = cluster.cat.categories.rename("cluster")
    return pd.Series(range(len(categories)), index=categories, name="x")


def get_whisker_location(ratios: pd.Series) -> pd.DataFrame:
    """Get the Y location for plotting p-value string.

    Y locations are the tip of the whisker.

    """
    _ratios = ratios.dropna()

    if len(_ratios) == 0:
        return 0

    lower, upper = np.quantile(_ratios, [0.25, 0.75])
    iqr = upper - lower
    upper_whisker = upper + 1.5 * iqr
    return upper_whisker


def munge_pvalues(input_file: str, x: pd.Series, y: pd.DataFrame) -> pd.DataFrame:
    """Make p-value table with locations for plotting."""
    return (
        pd.read_feather(input_file)
        .melt(id_vars="cluster", var_name="ratio_type", value_name="pvalue")
        .join(x, on="cluster")
        .merge(y, on=["cluster", "ratio_type"])
    )


if __name__ == "__main__":
    SNAKE = dict(
        ratios=snakemake.input.ratios,
        pvalues=snakemake.input.pvalues,
        output_prefix=snakemake.params[0],
    )

    # Debug Settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), "x-to-a-wf/scripts"))
    #     print(os.getcwd())
    # except:
    #     pass
    # snake = dict(
    #     ratios="../../output/x-to-a-wf/autosome_ratios_commonly_expressed_by_cell.feather",
    #     pvalues="../../output/x-to-a-wf/permuted_autosome_ratio_pvalues_commonly_expressed.feather",
    #     output_prefix=""
    # )

    main(SNAKE)
