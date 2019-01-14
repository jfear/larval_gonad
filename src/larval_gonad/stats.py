"""Module with helper functions for running statistical analysis.

This module houses various functions for helping run different types of statistical analysis.

"""
from typing import Tuple, Callable, Union

import numpy as np
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

