"""Combine coverage counts for each species.

Counts table is `YOgn x sample`

"""
from pathlib import Path

import pandas as pd


def main():
    df = pd.concat((read_data(file_name) for file_name in snakemake.input), axis=1, sort=True)
    df.reset_index().to_feather(snakemake.output[0])


def read_data(file_name):
    samplename = get_samplename(file_name)
    return (
        pd.read_csv(file_name, sep="\t", header=None, usecols=[0, 3], index_col=0)
        .rename_axis("YOgn")
        .squeeze()
        .rename(samplename)
    )


def get_samplename(file_name):
    return Path(file_name).stem


def debug():
    import os
    from larval_gonad.config import read_config
    from larval_gonad.mock import MockSnake

    try:
        os.chdir(os.path.join(os.getcwd(), "expression-atlas-wf"))
        print(os.getcwd())
    except:
        pass

    config = read_config("../config/common.yaml")
    tag = config["tag"]
    snakemake = MockSnake(
        input=[
            "../output/expression-atlas-wf/raw_counts/dana_AC_f_r1.tsv",
            "../output/expression-atlas-wf/raw_counts/dana_AC_m_r1.tsv",
            "../output/expression-atlas-wf/raw_counts/dana_DG_f_r1.tsv",
            "../output/expression-atlas-wf/raw_counts/dana_DG_f_r3.tsv",
        ]
    )
    return snakemake


# snakemake = debug()


if __name__ == "__main__":
    main()
