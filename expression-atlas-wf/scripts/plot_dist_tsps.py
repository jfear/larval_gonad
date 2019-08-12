"""Plot the distribution of TSPS.

TSPS is a measure of tissue specificity, where 0 is broad expression and 5 is
specific expression. I plot the distribution of TSPS in our data to determine
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
# INPUT_FILE = "../../output/expression-atlas-wf/dmel_tsps.feather"


def main():
    df = pd.read_feather(INPUT_FILE).set_index("FBgn")

    fig, ax = plt.subplots()
    sns.kdeplot(df.male_tsps.dropna(), label="male", ax=ax)
    sns.kdeplot(df.female_tsps.dropna(), label="female", ax=ax)
    ax.set(title="TSPS", ylabel="Density", xlabel="TSPS")
    plt.savefig(OUTPUT_FILE, bbox_inches="tight")


if __name__ == "__main__":
    main()
