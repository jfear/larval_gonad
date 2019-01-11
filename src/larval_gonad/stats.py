"""Module with helper functions for running statistical analysis.

This module houses various functions for helping run different types of statistical analysis.

"""
from typing import Tuple, Callable

import numpy as np
import scipy.stats


def permutation_sample(data1: np.ndarray, data2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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


