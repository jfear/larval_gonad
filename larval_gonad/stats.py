"""Module with helper functions for running statistical analysis.

This module houses various functions for helping run different types of statistical analysis.

"""
from typing import Tuple, Callable

import numpy as np
import scipy.stats


def permutation_sample(data1: np.ndarray, data2: np.ndarray) -> Tuple[np.ndarray]:
    """Generate a permuted sample.

    Combines two arrays, randomly mixes samples and split them back out into two arrays.

    Parameters
    ----------
    data1 :
        A numpy array of data.
    data2 :
        A numpy array of data.

    Returns
    -------
    tuple :
        A tuple of numpy.array objects that represent the permuted data1 and data2.

    """
    combined_data = np.concatenate((data1, data2))
    permuted_data = np.random.permutation(combined_data)
    permuted_sample_1 = permuted_data[:len(data1)]
    permuted_sample_2 = permuted_data[len(data1):]

    return permuted_sample_1, permuted_sample_2


def permuted_replicates(data1: np.ndarray, data2: np.ndarray, func: Callable, size: int = 1) -> np.ndarray:
    """Permute a number of samples and return a statistic.

    Takes two datasets, permutes them and returns the statistics from a function.

    Parameters
    ----------
    data1 : numpy.ndarray
        A numpy array of data.
    data2 : numpy.ndarray
        A numpy array of data.
    func : Callable
        A statistical function.
    size : int
        The number of permutations to perform.

    Returns
    -------
    numpy.ndarray :
        A numpy.array of p-values

    """
    replicates = np.empty(size)
    for i in range(size):
        permuted_sample_1, permuted_sample_2 = permutation_sample(data1, data2)
        replicates[i] = func(permuted_sample_1, permuted_sample_2)

    return replicates


def enrichment_statistic(data1: np.ndarray, data2: np.ndarray, test: str = 'ks'):
    """

    Parameters
    ----------
    data1 : numpy.ndarray
        A numpy array of data.
    data2 : numpy.ndarray
        A numpy array of data.
    test : str

    Returns
    -------
    float :
        The test statistic from `test`.

    """
    if test == 'ks':
        stat, _ = scipy.stats.ks_2samp(data1, data2)
    else:
        stat, _ = scipy.stats.mannwhitneyu(data1, data2)

    return stat
