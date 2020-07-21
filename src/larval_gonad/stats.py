"""Module with helper functions for running statistical analysis.

This module houses various functions for helping run different types of statistical analysis.

"""
from typing import List, Tuple, Union, Optional
from collections import namedtuple
from itertools import combinations
import concurrent.futures

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, norm, mannwhitneyu
from scipy.stats.contingency import margins
from statsmodels.stats.multitest import multipletests


def permute_sample(
    data1: np.ndarray, data2: np.ndarray, seed: Union[int, bool] = False
) -> Tuple[np.ndarray, np.ndarray]:
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
    permuted_sample_1 = permuted_data[: len(data1)]
    permuted_sample_2 = permuted_data[len(data1) :]

    return permuted_sample_1, permuted_sample_2


def calculate_permuted_pvalue(observed: float, permuted: list) -> float:
    """Calculate the permutation test p-value.

    Calculated as the number of permuted values more extreme (<=) the
    observed value.

    Args:
        observed (float): p-value from the statistical test.
        permuted (list): A list of p-values from the permutations.

    Returns:
        float: p-value as the proportion of more extreme observations.
    """
    n_more_extreme = np.sum(np.asarray(permuted) <= observed)
    return n_more_extreme / len(permuted)


PermutationResult = namedtuple(
    "PermutationResult", ["name1", "name2", "baseMean", "log2FoldChange", "p_value"]
)


class PairwisePermutationTest:
    """Run pairwise permutation tests."""

    def __init__(
        self,
        group_column: str,
        value_column: str,
        data: pd.DataFrame,
        method: Optional[str] = "mannwhitney",
        method_kws: Optional[dict] = None,
        order: Optional[List[str]] = None,
        n_permutations: int = 10_000,
        threads: int = 1,
    ):
        """Pairwise permutation tests.

        Args:
            group_column (str): The column name to group samples by.
            value_column (str): The column name with values.
            data (pd.DataFrame): A DataFrame containing `group_column`, `value_column`.
            method (str, optional): The statistical method to use for
            comparions. Must be one of ["mannwhitney", ].
            method_kws (dict, optional): Options to pass the method of choice.
            order (list, optional): The preferred order of the group columns.
            n_permutations (int): The number of permuations to perform. Defaults to 10,000
            threads (int, optional): Number of threads to run on. Defaults to 1.
        """
        self.group_column = group_column
        self.value_column = value_column
        self.data = data
        self.method = method
        self.method_kws = method_kws or {}
        self.n_permutations = n_permutations
        self.threads = threads

        self.ids = order or data[self.group_column].unique()
        self.comparisons = list(combinations(self.ids, 2))
        self._results = []

        self.mannwhitney_defaults = dict(use_continuity=True, alternative="two-sided")
        self.mannwhitney_defaults.update(self.method_kws)

    def fit(self):
        """Run the pairwise permutation test for all comparisons."""
        with concurrent.futures.ProcessPoolExecutor(self.threads) as executor:
            futures = []
            for name1, name2 in self.comparisons:
                data1 = self._get_data_from_group(name1)
                data2 = self._get_data_from_group(name2)

                if (len(data1) == 0) | (len(data2) == 0):
                    # No data, skip comparisons
                    continue

                futures.append(
                    executor.submit(self._run_comparison, name1, name2, data1, data2)
                )

            self._results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]

        return self

    def _get_data_from_group(self, name: str) -> np.ndarray:
        """Pull out rows corresponding to the group name."""
        return self.data.query(f"{self.group_column} == @name")[
            self.value_column
        ].values

    def _run_comparison(
        self, name1: str, name2: str, data1: np.ndarray, data2: np.ndarray
    ) -> PermutationResult:
        """Compare two distributions and calculate the p-value using permutation testing."""
        base_mean = np.mean(np.concatenate((data1, data2)))
        lfc = np.log2(np.mean(data1) / np.mean(data2))
        observed = self._run_test(data1, data2)
        permuted = self._permute_and_compare(data1, data2)
        p_value = calculate_permuted_pvalue(observed, permuted)
        return PermutationResult(name1, name2, base_mean, lfc, p_value)

    def _permute_and_compare(self, data1: np.ndarray, data2: np.ndarray) -> List[float]:
        """Run the actual permutation using multiple threads."""
        return [self._run_permutation(data1, data2) for _ in range(self.n_permutations)]

    def _run_test(self, *args) -> float:
        """Run the selected test statistic."""
        if self.method == "mannwhitney":
            return self._run_mannwhitney(*args)

    def _run_mannwhitney(self, data1: np.ndarray, data2: np.ndarray) -> float:
        """Wrapper around `scipy.stats.mannwhitneyu`.

        Args:
            data1 (np.ndarray): Passed as x to `scipy.stats.mannwhitneyu`.
            data2 (np.ndarray): Passed as y to `scipy.stats.mannwhitneyu`.

        Returns:
            float: The p-value from `scipy.stats.mannwhitneyu`.
        """
        _, pval = mannwhitneyu(data1, data2, **self.mannwhitney_defaults)
        return pval

    def _run_permutation(self, data1: np.ndarray, data2: np.ndarray) -> float:
        """Permutes the sample and runs the test.

        Returns:
            float: Value from `self._run_test` using permuted data.
        """
        return self._run_test(*permute_sample(data1, data2))

    @property
    def results(self):
        return (
            pd.DataFrame(self._results)
            .assign(
                name1=lambda x: pd.Categorical(
                    x.name1, categories=self.ids, ordered=True
                )
            )
            .assign(
                name2=lambda x: pd.Categorical(
                    x.name2, categories=self.ids, ordered=True
                )
            )
            .assign(padj_value=lambda x: multipletests(x.p_value, method="fdr_bh")[1])
        )


def adjusted_residuals(observed, expected):
    n = observed.sum().sum()
    rsum, csum = margins(observed)
    v = csum * rsum * (n - rsum) * (n - csum) / n ** 3
    return (observed - expected) / np.sqrt(v)


def make_big_table(
    obs, expected, resid, adj_resid, cell_chisqs, cell_qvals, cell_flags
):
    expected = pd.DataFrame(expected, index=obs.index, columns=obs.columns)
    cell_chisqs = pd.DataFrame(cell_chisqs, index=obs.index, columns=obs.columns)
    cell_qvals = pd.DataFrame(
        cell_qvals.reshape(resid.shape),
        index=adj_resid.index,
        columns=adj_resid.columns,
    )
    cell_flags = pd.DataFrame(
        cell_flags.reshape(resid.shape),
        index=adj_resid.index,
        columns=adj_resid.columns,
    )

    obs["type"] = "observed"
    obs = obs.set_index("type", append=True)

    expected["type"] = "expected"
    expected = expected.set_index("type", append=True)

    resid["type"] = "residual"
    resid = resid.set_index("type", append=True)

    adj_resid["type"] = "adj std residual"
    adj_resid = adj_resid.set_index("type", append=True)

    cell_chisqs["type"] = "X^2"
    cell_chisqs = cell_chisqs.set_index("type", append=True)

    cell_qvals["type"] = "fdr q-value"
    cell_qvals = cell_qvals.set_index("type", append=True)

    cell_flags["type"] = "flag_sig"
    cell_flags = cell_flags.set_index("type", append=True)

    _df = pd.concat(
        [obs, expected, resid, adj_resid, cell_chisqs, cell_qvals, cell_flags]
    ).reset_index(level="type")
    _df["type"] = pd.Categorical(
        _df["type"],
        ordered=True,
        categories=[
            "observed",
            "expected",
            "residual",
            "adj std residual",
            "X^2",
            "fdr q-value",
            "flag_sig",
        ],
    )
    return _df.set_index("type", append=True).sort_index().round(4)


def run_chisq(df, **kwargs):
    """ A helper function to run post hoc tests on chi^2.

    Parameters
    ----------
    df: a cross tabulation table
    kwargs: kwargs passed to statsmodels.stats.multitest.multipletests

    Returns
    -------
    pandas.DataFrame with chi-square and post-hoc tests.

    """
    obs = df.copy()
    stat, pval, degrees, expected = chi2_contingency(obs)
    print(f"𝛘^2: {stat:,.4f}, p-value: {pval:,.4f}, df: {degrees:,}")
    resid = obs - expected
    adj_resid = adjusted_residuals(df, expected)
    cell_chisqs = resid ** 2 / expected
    cell_flags, cell_qvals, _, _ = multipletests(
        norm.pdf(adj_resid).flatten(), method="fdr_bh", **kwargs
    )
    return make_big_table(
        obs, expected, resid, adj_resid, cell_chisqs, cell_qvals, cell_flags
    )

