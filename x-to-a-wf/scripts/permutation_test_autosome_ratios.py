import pandas as pd
from scipy.stats import mannwhitneyu

FNAME = snakemake.input[0]
ONAME = snakemake.output[0]

# Debug Settings
# FNAME = "output/x-to-a-wf/autosome_ratios_by_cell.feather"


def main():
    ratios = pd.read_feather(FNAME).set_index("cell_id")
    _ratios = ratios.copy()

    results = []
    for i in range(10_000):
        _ratios["cluster"] = _ratios.cluster.sample(frac=1).values
        for cluster, dd in _ratios.groupby("cluster"):
            obs = ratios.query(f'cluster == "{cluster}"')
            _, pval_x = mannwhitneyu(
                obs.x_to_a_ratio, dd.x_to_a_ratio, alternative="less"
            )
            _, pval_4 = mannwhitneyu(
                obs.fourth_to_a_ratio, dd.fourth_to_a_ratio, alternative="less"
            )
            _, pval_y = mannwhitneyu(
                obs.y_to_a_ratio, dd.y_to_a_ratio, alternative="greater"
            )
            results.append((cluster, pval_x <= 0.05, pval_4 <= 0.05, pval_y <= 0.05))

    df = pd.DataFrame(results, columns=["cluster", "pval_x", "pval_4", "pval_y"])
    pvals = 1 - df.groupby("cluster").mean().reindex(ratios.cluster.cat.categories).rename_axis('cluster')
    pvals.reset_index().to_feather(ONAME)


if __name__ == "__main__":
    main()

pd.Series.cat