from pickle import load
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

GONIA_V_CYTES = snakemake.input.gonia_vs_cytes
BACKGROUND_GENES = snakemake.input.background
INTRONLESS_GENES = snakemake.input.intronless

ONAME = snakemake.output[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
#     print(os.getcwd())
# except:
#     pass
# GONIA_V_CYTES = '../../output/seurat3-cluster-wf/combined_n3_gonia_vs_cytes.feather'
# BACKGROUND_GENES = '../../output/science_submission/background_fbgns.pkl'
# INTRONLESS_GENES = '../../output/science_submission/intron_less_genes.pkl'


def main():
    background = load(open(BACKGROUND_GENES, "rb"))
    intronless = load(open(INTRONLESS_GENES, "rb"))

    pct_expression = (
        pd.read_feather(GONIA_V_CYTES)
        .assign(pct_gonia=lambda x: x["pct.1"])
        .assign(pct_cyte=lambda x: x["pct.2"])
        .set_index("FBgn")
        .loc[:, ["pct_gonia", "pct_cyte"]]
        .reindex(background)
        .dropna()
    )
    pct_expression["intronless"] = np.where(
        pct_expression.index.isin(intronless), "intronless", "has_intron"
    )

    plt.style.use("scripts/figure_styles.mplstyle")
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=plt.figaspect(1 / 2))
    sns.boxplot(
        "intronless",
        "pct_gonia",
        data=pct_expression,
        ax=ax1,
        linewidth=0.5,
        showfliers=False,
        notch=True,
    )
    ax1.set(title="SP Expression", xlabel="", ylabel="Percent of Cells", ylim=(-0.1, 1.1))

    sns.boxplot(
        "intronless",
        "pct_cyte",
        data=pct_expression,
        ax=ax2,
        linewidth=0.5,
        showfliers=False,
        notch=True,
    )
    ax2.set(title="Spermatocyte Expression", xlabel="", ylabel="")

    fig.savefig(ONAME, bbox_inches="tight")


if __name__ == "__main__":
    main()
