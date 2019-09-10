import os
import re
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import cellranger_umi


def main():
    REP = re.findall(r"testis(\d)", snakemake.wildcards.sample)[0]
    df = get_data(REP)

    fig, ax = plt.subplots()
    defaults = dict(x="cell_count", y="UMI", kind="scatter", rasterized=True, ax=ax)
    df.plot(label="Empty", color="lightgray", **defaults)
    df[df.is_cell].plot(label="Cell", color="k", **defaults)
    df[df.is_doublet].plot(label="Doublet", color="r", marker="^", **defaults)
    ax.set(yscale="log", xscale="log", title=f"Testis {REP}")
    ax.legend()
    sns.despine(ax=ax)

    fig.savefig(snakemake.output[0])


def get_data(rep):
    cell_calls = pd.read_feather(
        snakemake.input.cell_calls, columns=["cell_id", "is_cell"]
    ).set_index("cell_id")

    doublet_calls = [x.strip() for x in open(snakemake.input.doublets)]
    cell_calls["is_doublet"] = np.where(cell_calls.index.isin(doublet_calls), True, False)

    umi = cellranger_umi(snakemake.input.umi).groupby("cell_id").size().rename("UMI")
    umi.index = f"rep{rep}_" + umi.index

    cell_calls = cell_calls.reindex(umi.index).fillna(False)
    cell_calls["UMI"] = umi.values

    cell_calls_sorted = cell_calls.sort_values("UMI", ascending=False)
    cell_calls_sorted["cell_count"] = np.arange(1, cell_calls_sorted.shape[0] + 1)

    return cell_calls_sorted


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        SAMPLE = "testis1"
        snakemake = snakemake_debug(
            workdir="cellselection-wf",
            input=dict(
                cell_calls=f"../output/cellselection-wf/{SAMPLE}_combined_cell_calls.feather",
                doublets=f"../output/cellselection-wf/{SAMPLE}_scrublet_dublets.txt",
                umi=f"../output/cellranger3-wf/{SAMPLE}/outs/molecule_info.h5",
            ),
            wildcards=dict(sample=SAMPLE),
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
