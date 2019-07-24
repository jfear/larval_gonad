import re
from pathlib import Path

import pandas as pd

COUNT_FILES = snakemake.input.counts
ORTHROLOG_FILES = snakemake.input.orthologs
OUTPUT_FILE = snakemake.output[0]

# Debug settings
# import os
# os.chdir('neox-wf/scripts')
# COUNT_FILES = [
#     '../../output/neox-wf/raw_counts/dana_AC_f_r1.tsv',
#     '../../output/neox-wf/raw_counts/dana_AC_f_r2.tsv',
#     '../../output/neox-wf/raw_counts/dana_AC_f_r3.tsv',
#     '../../output/neox-wf/raw_counts/dper_AC_f_r1.tsv',
#     '../../output/neox-wf/raw_counts/dper_AC_f_r2.tsv',
#     '../../output/neox-wf/raw_counts/dper_AC_f_r3.tsv'
# ]
# ORTHROLOG_FILES = [
#     '../../output/neox-wf/orthologs/dana.tsv',
#     '../../output/neox-wf/orthologs/dmel.tsv',
#     '../../output/neox-wf/orthologs/dmoj.tsv',
#     '../../output/neox-wf/orthologs/dper.tsv',
#     '../../output/neox-wf/orthologs/dpse.tsv',
#     '../../output/neox-wf/orthologs/dvir.tsv',
#     '../../output/neox-wf/orthologs/dwil.tsv',
#     '../../output/neox-wf/orthologs/dyak.tsv'
# ]


def main():
    orthologs = (
        pd.concat((read_orthologs(fname) for fname in ORTHROLOG_FILES), sort=True)
        .unstack(level="species")
        .rename_axis("FBgn")
        .reset_index()
    )

    counts = pd.concat(
        (read_counts(fname, orthologs) for fname in COUNT_FILES), sort=True, axis=1
    ).rename_axis("FBgn")

    counts.reset_index().to_feather(OUTPUT_FILE)


def read_orthologs(fname):
    return (
        pd.read_csv(fname, sep="\t", index_col="Dmel", na_values="-", usecols=["Dmel", "YOgnID"])
        .assign(species=re.findall(r"(d\w+)\.tsv", fname)[0])
        .set_index("species", append=True)
        .squeeze()
    )


def read_counts(fname, orthologs):
    name = Path(fname).stem
    species = name.split("_")[0]

    if (species == "orgR") or (species == "w1118"):
        species = "dmel"

    df = pd.read_csv(
        fname,
        sep="\t",
        header=None,
        names=[species, "jaccard", "FBgn", "cnt"],
        na_values="-",
        usecols=[species, "cnt"],
    )

    return df.merge(orthologs[["FBgn", species]], on=species).set_index("FBgn").cnt.rename(name)


if __name__ == "__main__":
    main()
