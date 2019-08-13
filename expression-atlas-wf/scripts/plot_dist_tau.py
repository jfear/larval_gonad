"""Plot the distribution of tau.

Tau is a measure of tissue specificity, where 0 is broad expression and 1 is
specific expression. I plot the distribution of Tau in our data to determine
a valid cutoff.

"""
import matplotlib

matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

INPUT_FILE = snakemake.input[0]
OUTPUT_FILE = snakemake.output[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'expression-atlas-wf/scripts'))
#     print(os.getcwd())
# except:
#     pass
# INPUT_FILE = "../../output/expression-atlas-wf/dmel_tau.feather"


def main():
    df = pd.read_feather(INPUT_FILE).set_index("FBgn")

    fig, ax = plt.subplots()
    sns.kdeplot(df.male_tau.dropna(), label="male", ax=ax)
    sns.kdeplot(df.female_tau.dropna(), label="female", ax=ax)
    ax.set(title="Tau", ylabel="Density", xlabel="tau")
    plt.savefig(OUTPUT_FILE, bbox_inches="tight")


if __name__ == "__main__":
    main()
