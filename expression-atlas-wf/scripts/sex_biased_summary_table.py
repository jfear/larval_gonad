import re

import pandas as pd

from larval_gonad.io import pickle_load


def main():
    res = [munge(file_name) for file_name in SNAKE["input"]]
    df = pd.DataFrame(res, columns=["direction", "species", "tissue", "count"])
    df.to_feather(SNAKE["output"])


def munge(file_name):
    attrs = re.findall(r".*/sex_biased_(\w+)_(\w+)_(\w+)_fbgns.pkl", file_name)[0]
    cnt = len(pickle_load(file_name))
    return [*attrs, cnt]


if __name__ == "__main__":
    SNAKE = dict(input=snakemake.input, output=snakemake.output[0])

    # Debug settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), "expression-atlas-wf/scripts"))
    #     print(os.getcwd())
    # except:
    #     pass

    # SNAKE = dict(
    #     input=["../../output/expression-atlas-wf/sex_biased_male_orgR_GO_fbgns.pkl"], output=""
    # )

    main()
