from pathlib import Path
from itertools import product
from pathlib import Path

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import metrics

# Globals
SAMPLE = snakemake.wildcards["sample"]
INPUTS = snakemake.input
OUTPUT = snakemake.output[0]


def main():
    # Get all cluster calls
    df = pd.concat((get_data(x) for x in INPUTS), axis=1, sort=True).dropna()

    # Calculate the pairwise adjusted RAND Score
    res = (
        pd.DataFrame(
            [
                (c1, c2, metrics.adjusted_rand_score(df[c1], df[c2]))
                for c1, c2 in product(df.columns, df.columns)
            ],
            columns=["c1", "c2", "AdjRand"],
        )
        .set_index(["c1", "c2"])
        .sort_index()
        .unstack()
    )
    res.columns = res.columns.droplevel(0)

    # Get row and col ordering based on numeric representation of names
    idx_order = (
        pd.DataFrame(index=res.index)
        .assign(low=lambda x: x.index.str.extract("^(\d+)_.*", expand=False))
        .assign(high=lambda x: x.index.str.extract("^.*?_(\d+|None)_.*", expand=False))
        .assign(
            mito=lambda x: x.index.str.extract("^.*?_.*?_(\d+|None)_.*", expand=False)
        )
        .assign(
            ribo=lambda x: x.index.str.extract("^.*?_.*?_.*?_(\d+|None)", expand=False)
        )
        .applymap(make_numeric)
        .sort_values(["low", "high", "mito", "ribo"])
        .index
    )

    # Plot cluster map of RAND Scores
    fig = plt.figure(figsize=(20, 20))
    ax = sns.heatmap(
        res.reindex(idx_order).reindex(idx_order.values, axis=1),
        xticklabels=True,
        yticklabels=True,
        square=True,
        vmin=0.5,
        vmax=1,
        cbar_kws={"label": "Adj Rand Score"},
    )
    ax.set(
        xlabel="Clusters (low gene _ high gene _ mito _ ribo)",
        ylabel="Clusters (low gene _ high gene _ mito _ ribo)",
        title=SAMPLE,
    )
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

    plt.savefig(OUTPUT, bbox_inches="tight")


def get_data(fname):
    name = "_".join(Path(fname).stem.split("_")[1:5])

    return (
        pd.read_feather(fname, columns=["cell_id", "cluster"])
        .rename(columns={"cluster": name})
        .set_index("cell_id")
        .squeeze()
    )


def make_numeric(x):
    if x == "None":
        return np.nan
    else:
        return int(x)


if __name__ == "__main__":
    main()
