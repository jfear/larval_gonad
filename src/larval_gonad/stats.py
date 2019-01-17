"""Module with helper functions for running statistical analysis.

This module houses various functions for helping run different types of statistical analysis.

"""
from typing import Tuple, Callable, Union
from collections import namedtuple

import numpy as np
import pandas as pd
import scipy.stats


def permutation_sample(data1: np.ndarray, data2: np.ndarray, seed: Union[int, bool] = False) \
        -> Tuple[np.ndarray, np.ndarray]:
    """Generate a permuted sample.

    Combines two arrays, randomly mixes samples and split them back out into two arrays.

    Parameters
    ----------
    data1 :
        A numpy array of data.
    data2 :
        A numpy array of data.
    seed : int
        Set a random seed.

    Returns
    -------
    tuple :
        A tuple of numpy.array objects that represent the permuted data1 and data2.

    """
    if seed:
        np.random.seed(seed)

    combined_data = np.concatenate((data1, data2))
    permuted_data = np.random.permutation(combined_data)
    permuted_sample_1 = permuted_data[:len(data1)]
    permuted_sample_2 = permuted_data[len(data1):]

    return permuted_sample_1, permuted_sample_2


def _drop_zeros(chrom1, chrom2):
    """Helper function to remove pairs of indexes where chrom2 is zero.

    If chrom2 contains zeros then we end up with a zero division error. The easiest solution is to drop indices from
    chrom1 and chrom2 where chrom2 == 0. This function returns the two arrays back where (chrom2 != 0).

    """
    chrom1 = np.array(chrom1)
    chrom2 = np.array(chrom2)
    non_zero_ids = np.where(chrom2 > 0)
    return chrom1[non_zero_ids], chrom2[non_zero_ids]


def permutation_test_chrom1_lt_chrom2(target_chrom: np.ndarray, autosome: np.ndarray, size: int = 1_000, ) -> float:
    """Calculates if the median chromosome ratio is extreme.

    Calculates the median chromosomal ratio and compares to a permutation set. A p-value <=0.05 indicates that
    chrom1 < chrom2.

    Parameters
    ----------
    target_chrom : np.ndarray
        An array-like list or read counts from chromosome 1.
    autosome : np.ndarray
        An array-like list or read counts from chromosome 2.
    size : int
        The number of permutations to run.

    Returns
    -------
    float
        The probability of having a more extreme median ratio.

    """
    target_chrom, autosome = _drop_zeros(target_chrom, autosome)

    func = lambda c1, c2: np.median(c1 / c2)
    observed_median_ratio = func(target_chrom, autosome)
    permutation_results = np.empty(size)
    for i in range(size):
        permuted_chrom1, permuted_chrom2 = _drop_zeros(*permutation_sample(target_chrom, autosome))
        permutation_results[i] = func(permuted_chrom1, permuted_chrom2)
    p_value = sum(permutation_results <= observed_median_ratio) / len(permutation_results)
    return p_value


MannWhitneyResult = namedtuple('MannWhitneyResult', 'cell_id flag_X_depleted')


def mannwhitneyu_cell_level_x_to_autosome(chromosome: str, umi: str, data: pd.DataFrame, cell_id: str,
                                          p_value_cutoff: float = 0.05) -> MannWhitneyResult:
    """Calculates the Mann-Whitney U statsitic comparing median X vs median Autosome expression.

    Take data from a single-cell and compares median X vs median A of expressed genes.

    Parameters
    ----------
    chromosome : str
        The name of the column that contains chromosome labels. Labels must be `X` or `A`.
    umi : str
        The name of the column that contains UMI counts.
    data : pandas.DataFrame
        Data from a single-cell indexed by gene_id with UMI counts and chromosome annotations.
    cell_id : str
        The cell id being analyzed, will be included in the return results.
    p_value_cutoff : float
        P-value cutoff to use.

    Returns
    -------
    MannWhitneyResult
        A namedtuple with (cell_id, flag_X_depleted)

    """
    if data[chromosome].unique().shape[0] == 1:
        return MannWhitneyResult(cell_id, np.nan)

    x_genes = data.query(f'{chromosome} == "X"')[umi].value
    a_genes = data.query(f'{chromosome} == "A"')[umi].value

    if x_genes.shape[0] < 100 and a_genes.shape[0] < 100:
        return MannWhitneyResult(cell_id, np.nan)

    stat, p_value = scipy.stats.mannwhitneyu(x_genes, a_genes, alternative='less')

    if p_value < p_value_cutoff:
        return MannWhitneyResult(cell_id, True)

    return MannWhitneyResult(cell_id, False)
