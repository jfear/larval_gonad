import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import shelve_load
from larval_gonad.plotting import format_pval


def main():
    db = shelve_load(SNAKE["input_file"], "data", "pvalues")
    df = db["data"]
    pvals = db["pvalues"]

    ax = sns.barplot(
        "chrom", "Proportion of Genes", data=df, order=SNAKE["chrom_order"][:-1], color="C0"
    )
    pvals.query("chrom != 'Y'").apply(lambda x: format_pval(ax, x.x, x.y, x.qval), axis=1)
    ax.set(ylim=(0, 1), title=SNAKE['title'])
    
    plt.savefig(SNAKE['output_file'])


if __name__ == "__main__":
    SNAKE = dict(
        input_file=snakemake.input[0],
        output_file=snakemake.output[0],
        chrom_order=snakemake.params.chrom_order,
        title=snakemake.params.title,
    )
    plt.style.use("scripts/figure_styles.mplstyle")

    # Debug Settings
    # import os

    # try:
    #     os.chdir(os.path.join(os.getcwd(), "science_submission/scripts"))
    #     print(os.getcwd())
    # except:
    #     pass
    # from larval_gonad.config import read_config

    # config = read_config("../../config/common.yaml")
    # SNAKE = dict(
    #     input_file="../../output/science_submission/db/demasculinization.larval_testis.dat",
    #     output_file="",
    #     chrom_order=config["chrom_order"],
    #     title="Larval Testis-Biased Gene Expression",
    # )
    # plt.style.use("figure_styles.mplstyle")

    main()
