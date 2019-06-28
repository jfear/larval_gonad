"""Library for normalization of RNA-seq data."""
import numpy as np


def cpm(df, scale=1e6, log=None):
    """Coverage per million.

    Normalizes data by dividing by the library size and multiplying by scaling
    factor (typically 1 million).

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with genes as rows and samples as columns.
    scale : int or float
        Scaling factor to multiple counts.
    log : None or function or str
        If a function is giving will apply this to the data before returning.
        If a string is given then it uses numpy log functions that correspond
        to the provided string.

    """
    totals = df.sum()
    if log is None:
        log = lambda x: x - 1
    elif log == "log2":
        log = np.log2
    elif log == "log10":
        log = np.log10
    elif log == "ln":
        log = np.log

    return log(df / (totals / scale) + 1)


def rpkm(df, gene_length, scale_library=1e6, scale_length=1e3, log=None):
    """Reads per Million mapped reads per kilobase gene mode.

    Calcualtes RPKM which normalizes by library size and gene length.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with genes as rows and samples as columns.
    gene_length : pd.Series
        Series were index matches df.index and values are gene lengths.
    scale_library : int or float
        Scaling factor to scale the library size.
    scale_length : int or float
        Scaling factor to scale gene model lengths.
    log : None or function or str
        If a function is giving will apply this to the data before returning.
        If a string is given then it uses numpy log functions that correspond
        to the provided string.

    """
    totals = df.sum()
    if log is None:
        log = lambda x: x - 1
    elif log == "log2":
        log = np.log2
    elif log == "log10":
        log = np.log10
    elif log == "ln":
        log = np.log

    return log((df / (totals / scale_library)).T / (gene_length / scale_length) + 1).T


def tpm(df, gene_length, scale_library=1e6, scale_length=1e3, log=None):
    """Transcripts Per Killobase Million.

    Calcualtes TPM which normalizes by library size and gene length.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with genes as rows and samples as columns.
    gene_length : pd.Series
        Series were index matches df.index and values are gene lengths.
    scale_library : int or float
        Scaling factor to scale the library size.
    scale_length : int or float
        Scaling factor to scale gene model lengths.
    log : None or function or str
        If a function is giving will apply this to the data before returning.
        If a string is given then it uses numpy log functions that correspond
        to the provided string.

    """

    rpk = (df.T / (gene_length / scale_length)).T
    totals = rpk.sum()

    if log is None:
        log = lambda x: x - 1
    elif log == "log2":
        log = np.log2
    elif log == "log10":
        log = np.log10
    elif log == "ln":
        log = np.log

    return log((rpk / (totals / scale_library)) + 1)


def zscore(df):
    mu = df.mean(axis=1)
    sigma = df.std(axis=1)
    return df.subtract(mu, axis="rows").div(sigma, axis="rows")
