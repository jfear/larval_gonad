import pytest
import numpy as np

from larval_gonad.stats import permutation_sample, enrichment_statistic, permuted_replicates
np.random.seed(42)


@pytest.fixture(scope='session')
def example_data():
    data1 = np.array([1, 2, 3, 4, 5])
    data2 = np.array([11, 12, 13, 14])
    return data1, data2


def test_permutation_sample(example_data):
    data1, data2 = example_data
    permuted_sample_1, permuted_sample_2 = permutation_sample(data1, data2)
    assert len(permuted_sample_1) == len(data1)
    assert np.all(permuted_sample_1 == np.array([13, 2, 11, 1, 14]))
    assert len(permuted_sample_2) == len(data2)
    assert np.all(permuted_sample_2 == np.array([3, 5, 4, 12]))


def test_enrichment_statistic_ks():
    stat1 = enrichment_statistic(
        np.array([1, 2, 3, 4, 5]),
        np.array([100, 110, 120, 130, 140])
    )

    stat2 = enrichment_statistic(
        np.array([1, 2, 3, 4, 5]),
        np.array([2, 3, 4, 5, 6])
    )

    assert stat1 > stat2


def test_enrichment_statistic_mannwhiteneyu():
    stat1 = enrichment_statistic(
        np.array([1, 2, 3, 4, 5]),
        np.array([100, 110, 120, 130, 140]),
        'mannwhiteneyu'
    )

    stat2 = enrichment_statistic(
        np.array([1, 2, 3, 4, 5]),
        np.array([2, 3, 4, 5, 6])
    )

    assert stat1 < stat2


def test_permuted_replicates(example_data):
    data1, data2 = example_data
    obs = enrichment_statistic(data1, data2)
    replicates = permuted_replicates(data1, data2, func=enrichment_statistic, size=5)
    p_value = np.sum(replicates > obs) / len(replicates)
    assert len(replicates) == 5
    assert p_value == 0


