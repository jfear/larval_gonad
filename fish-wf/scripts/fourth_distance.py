import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

from larval_gonad import plotting # pylint: disable=unused-import

plt.style.use("minimal")


def main():
    distance = (
        pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
        .query("Phase == 'middle/late'")
        .iloc[:, 6:8]
    )
    distance.columns = ["X", "4"]
    distance_melted = distance.melt(var_name="Chromosome", value_name="Distance (μm)")

    ax = sns.boxplot(
        "Chromosome",
        "Distance (μm)",
        data=distance_melted,
        palette=snakemake.params.colors,
        notch=True,
        showfliers=False,
    )

    # Clean up plot
    ax.set(xlabel="")
    sns.despine(ax=ax)

    # Test that not significant
    _, pval = wilcoxon(
        distance.X, distance["4"], alternative="two-sided", correction=True
    )
    if pval <= 0.05:
        # Extend axis and add NS.
        _max = distance.max().max() + .1
        ax.set_ylim(None, _max)

        ax.text(
            0.5, 0.99, f"p = {pval:.2E}", transform=ax.transAxes, va="top", ha="center"
        )
        l = plt.Line2D(
            [0.3, 0.7], [0.94, 0.94], transform=ax.transAxes, color="k", lw=0.8, ls="-"
        )
        ax.add_line(l)

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    main()
