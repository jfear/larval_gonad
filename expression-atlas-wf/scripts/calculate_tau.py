import numpy as np
import pandas as pd

from larval_gonad.normalization import tpm

COUNTS_FILE = snakemake.input["counts"]
SAMPLETABLE_FILE = snakemake.input["sampletable"]
GENE_ANNOT = snakemake.input["gene_annot"]
OUTPUT_FILE = snakemake.output[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'expression-atlas-wf/scripts'))
#     print(os.getcwd())
# except:
#     pass
# COUNTS_FILE = '../../output/expression-atlas-wf/raw_counts.feather'
# SAMPLETABLE_FILE = '../config/sampletable.tsv'
# GENE_ANNOT = '../../references/gene_annotation_dmel_r6-26.feather'


def main():
    df = read_counts_data()

    male_tpm = split_and_normalize(df, "male")
    female_tpm = split_and_normalize(df, "female")

    df_tau = pd.concat(
        [
            male_tpm.apply(tau, axis=1).rename("male_tau"),
            female_tpm.apply(tau, axis=1).rename("female_tau"),
        ],
        axis=1,
        sort=True,
    ).rename_axis("FBgn")

    df_tau.reset_index().to_feather(OUTPUT_FILE)


def read_counts_data():
    sampletable = (
        pd.read_csv(
            SAMPLETABLE_FILE, sep="\t", usecols=["samplename", "species", "tissue", "sex", "rep"]
        )
        .query("species == 'Drosophila melanogaster'")
        .set_index("samplename")
        .drop("species", axis=1)
    )

    df = (
        pd.read_feather(COUNTS_FILE, columns=["FBgn"] + sampletable.index.tolist())
        .set_index("FBgn")
        .dropna()
        .pipe(lambda x: x[x.sum(axis=1) > 0])
        .T.join(sampletable[["tissue", "sex"]])
        .groupby(["tissue", "sex"])
        .sum()
    )

    return df


def split_and_normalize(df, sex):
    fbgn2length = (
        pd.read_feather(GENE_ANNOT, columns=["FBgn", "length"]).set_index("FBgn").squeeze()
    )
    split_df = split_sex(df, sex)
    return tpm(split_df, fbgn2length).dropna()


def split_sex(df, sex):
    sexed_df = df.loc[(slice(None), sex), :].T
    sexed_df.columns = sexed_df.columns.droplevel(-1)
    return sexed_df


def tau(x: np.ndarray):
    """Calculate the tissue specificity measure tau as in:

    > Yanai, Itai, Hila Benjamin, Michael Shmoish, Vered Chalifa-Caspi, Maxim
    > Shklar, Ron Ophir, Arren Bar-Even, et al. 2005. “Genome-Wide Midrange
    > Transcription Profiles Reveal Expression Level Relationships in Human
    > Tissue Specification.” Bioinformatics 21 (5): 650–59.

    Example
    -------
    >>> tau(np.array([1, 1, 1]))
    0.0
    >>> tau(np.array([0, 0, 0]))
    np.nan
    >>> tau(np.array([1, 0, 0]))
    1.0

    """
    _max = x.max()
    if _max == 0:
        return np.nan
    _n = x.shape[0]
    xhat = x / _max
    tau = (1 - xhat).sum() / (_n - 1)
    return tau


if __name__ == "__main__":
    main()
