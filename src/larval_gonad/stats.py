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


def chromosome_permutation_test(chrom1: np.ndarray, chrom2: np.ndarray, size: int = 1_000, ) -> float:
    """Calculates if the median chromosome ratio is extreme.

    Calculates the median chromosomal ratio and compares to a permutation set. A p-value <=0.05 indicates that
    chrom1 < chrom2.

    Parameters
    ----------
    chrom1 : np.ndarray
        An array-like list or read counts from chromosome 1.
    chrom2 : np.ndarray
        An array-like list or read counts from chromosome 2.
    size : int
        The number of permutations to run.

    Returns
    -------
    float
        The probability of having a more extreme median ratio.

    """
    func = lambda c1, c2: np.median(c1 / c2)
    observed_median_ratio = func(chrom1, chrom2)
    permutation_results = np.empty(size)
    for i in range(size):
        permuted_chrom1, permuted_chrom2 = permutation_sample(chrom1, chrom2)
        permutation_results[i] = func(permuted_chrom1, permuted_chrom2)
    p_value = sum(permutation_results <= observed_median_ratio) / len(permutation_results)
    return p_value

