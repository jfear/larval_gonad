import re
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import cellranger_umi

CELL_CALLS = snakemake.input["cell_calls"]
DOUBLET_CALLS = snakemake.input["doublets"]
UMI_COUNTS = snakemake.input["umi"]
OUTPUT_FILE = snakemake.output[0]

SAMPLE = snakemake.wildcards.sample
REP = re.findall(r"testis(\d)", SAMPLE)[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
#     print(os.getcwd())
# except:
#     pass

# SAMPLE = 'testis1'
# REP = re.findall(r"testis(\d)", SAMPLE)[0]

# CELL_CALLS = f'../../output/cellselection-wf/{SAMPLE}_combined_cell_calls.feather'
# DOUBLET_CALLS = f'../../output/cellselection-wf/{SAMPLE}_scrublet_dublets.txt'
# UMI_COUNTS = f'../../output/cellranger3-wf/{SAMPLE}/outs/molecule_info.h5'


def main():
    df = get_data()

    plt.style.use("scripts/figure_styles.mplstyle")
    fig, ax = plt.subplots()
    defaults = dict(x="cell_count", y="UMI", kind="scatter", rasterized=True, ax=ax)
    df.plot(label="Empty", color="lightgray", **defaults)
    df[df.is_cell].plot(label="Cell", color="k", **defaults)
    df[df.is_doublet].plot(label="Doublet", color="r", marker="^", **defaults)
    ax.set(yscale="log", xscale="log", title=f"Testis {REP}")
    ax.legend()
    sns.despine(ax=ax)

    fig.savefig(OUTPUT_FILE, bbox_inches="tight")


def get_data():
    cell_calls = pd.read_feather(CELL_CALLS, columns=["cell_id", "is_cell"]).set_index("cell_id")

    doublet_calls = [x.strip() for x in open(DOUBLET_CALLS)]
    cell_calls["is_doublet"] = np.where(cell_calls.index.isin(doublet_calls), True, False)

    umi = cellranger_umi(UMI_COUNTS).groupby("cell_id").size().rename("UMI")
    umi.index = f"rep{REP}_" + umi.index

    cell_calls = cell_calls.reindex(umi.index).fillna(False)
    cell_calls["UMI"] = umi.values

    cell_calls_sorted = cell_calls.sort_values("UMI", ascending=False)
    cell_calls_sorted["cell_count"] = np.arange(1, cell_calls_sorted.shape[0] + 1)

    return cell_calls_sorted


if __name__ == "__main__":
    main()
