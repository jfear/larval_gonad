"""Boxplot of GSEA Enrichment Scores."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import larval_gonad.plotting

plt.style.use("minimal")

def main():
    df = pd.read_feather(snakemake.input[0])

    ax: plt.Axes
    ax = sns.boxplot(
        x="cluster",
        y="male",
        data=df,
        palette=snakemake.params.colors,
        order=snakemake.params.order,
        notch=True,
        showfliers=False,
    )
    sns.despine(ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("GSEA Enrichment Score")
    new_labels = [
        snakemake.params.names[label.get_text()]
        for label in ax.get_xticklabels()
    ]
    ax.set_xticklabels(new_labels, rotation=45)
    
    plt.savefig(snakemake.output[0])



if __name__ == "__main__":
    main()
