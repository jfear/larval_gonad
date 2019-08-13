"""Permute autosome ratios and calculate p-value.

Use a permutation test to calculate a p-values for comparing autosome ratios
among clusters.

1. Permute cluster ids 10,000 times
2. Compare distributions of autosome ratios between permuted and observed values.
3. Count the number of times observed values were more extreme than permuted values.
4. Calculate p-value based on 3.

"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def main(snake):
    ratios = pd.read_feather(snake["input_file"]).set_index("cell_id")
    _ratios = ratios.copy()

    results = []
    for i in range(10_000):
        _ratios["cluster"] = _ratios.cluster.sample(frac=1).values
        for cluster, dd in _ratios.groupby("cluster"):
            obs = ratios.query(f'cluster == "{cluster}"')

            perm_values_x = run_mannwhitney(obs.x_to_a_ratio, dd.x_to_a_ratio)
            perm_values_4 = run_mannwhitney(obs.fourth_to_a_ratio, dd.fourth_to_a_ratio)
            perm_values_y = run_mannwhitney(
                obs.y_to_a_ratio, dd.y_to_a_ratio, alternative="greater"
            )

            results.append(
                (
                    cluster,
                    summarize_permutation(perm_values_x),
                    summarize_permutation(perm_values_4),
                    summarize_permutation(perm_values_y),
                )
            )

    df = pd.DataFrame(
        results, columns=["cluster", "x_to_a_ratio", "fourth_to_a_ratio", "y_to_a_ratio"]
    )

    pvalue = 1 - (
        df.groupby("cluster")
        .apply(lambda x: x.mean())
        .reindex(ratios.cluster.cat.categories)
        .rename_axis("cluster")
    )

    pvalue.reset_index().to_feather(snake["output_file"])


def run_mannwhitney(obs: np.ndarray, permuted: np.ndarray, alternative="less") -> int:
    """Compare distributions of observed and permuted values."""
    _obs = obs.dropna()
    _permuted = permuted.dropna()
    if len(_obs) == 0:
        return np.nan

    _, pval = mannwhitneyu(_obs, _permuted, alternative=alternative)
    return pval


def summarize_permutation(perm_value: int, alpha=0.05) -> bool:
    """Determine the number of times the observed was more extreme."""
    if np.isnan(perm_value):
        return np.nan
    return perm_value <= alpha


if __name__ == "__main__":
    SNAKE = dict(input_file=snakemake.input[0], output_file=snakemake.output[0])

    # Debug Settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), 'x-to-a-wf/scripts'))
    #     print(os.getcwd())
    # except:
    #     pass
    # SNAKE = dict(
    #     input_file="../../output/x-to-a-wf/autosome_ratios_commonly_expressed_by_cell.feather",
    #     output_file=''
    # )

    main(SNAKE)
