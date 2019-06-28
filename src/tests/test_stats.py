import pytest
import numpy as np

from larval_gonad.stats import permutation_sample, permutation_test_chrom1_lt_chrom2

np.random.seed(42)

SAMPLE_SIZE = 1_000


@pytest.fixture(scope="session")
def toy_data():
    data1 = np.array([1, 2, 3, 4, 5])
    data2 = np.array([11, 12, 13, 14])
    return data1, data2


@pytest.fixture(scope="session")
def random_counts():
    return np.random.lognormal(6, size=SAMPLE_SIZE)


@pytest.fixture(scope="session")
def chrom1_chrom2_large_difference(random_counts):
    chrom1 = random_counts
    chrom2 = chrom1 + np.random.normal(500, size=SAMPLE_SIZE)
    return chrom1, chrom2


@pytest.fixture(scope="session")
def chrom1_chrom2_small_difference(random_counts):
    chrom1 = random_counts
    chrom2 = chrom1 + np.random.normal(50, size=SAMPLE_SIZE)
    return chrom1, chrom2


@pytest.fixture(scope="session")
def chrom1_chrom2_minute_difference(random_counts):
    chrom1 = random_counts
    chrom2 = chrom1 + np.random.normal(5, size=SAMPLE_SIZE)
    return chrom1, chrom2


@pytest.fixture(scope="session")
def chrom1_chrom2_no_difference(random_counts):
    chrom1 = random_counts
    chrom2 = chrom1 + np.random.normal(0, size=SAMPLE_SIZE)
    return chrom1, chrom2


def test_permutation_sample(toy_data):
    data1, data2 = toy_data
    permuted_sample_1, permuted_sample_2 = permutation_sample(data1, data2, seed=42)
    assert len(permuted_sample_1) == len(data1)
    assert np.all(permuted_sample_1 == np.array([13, 2, 11, 1, 14]))
    assert len(permuted_sample_2) == len(data2)
    assert np.all(permuted_sample_2 == np.array([3, 5, 4, 12]))


def test_chrom1_lt_chrom2_large_difference(chrom1_chrom2_large_difference):
    chrom1, chrom2 = chrom1_chrom2_large_difference
    assert permutation_test_chrom1_lt_chrom2(chrom1, chrom2, size=100) <= 0.05


def test_chrom1_lt_chrom2_small_difference(chrom1_chrom2_small_difference):
    chrom1, chrom2 = chrom1_chrom2_small_difference
    assert permutation_test_chrom1_lt_chrom2(chrom1, chrom2, size=100) <= 0.05


def test_chrom1_lt_chrom2_minute_difference(chrom1_chrom2_minute_difference):
    chrom1, chrom2 = chrom1_chrom2_minute_difference
    assert permutation_test_chrom1_lt_chrom2(chrom1, chrom2, size=100) > 0.05


def test_chrom1_lt_chrom2_no_difference(chrom1_chrom2_no_difference):
    chrom1, chrom2 = chrom1_chrom2_no_difference
    assert permutation_test_chrom1_lt_chrom2(chrom1, chrom2, size=100) > 0.05


def test_chrom1_gt_chrom2_large_difference(chrom1_chrom2_large_difference):
    chrom1, chrom2 = chrom1_chrom2_large_difference
    assert permutation_test_chrom1_lt_chrom2(chrom2, chrom1, size=100) > 0.05


def test_chrom1_gt_chrom2_small_difference(chrom1_chrom2_small_difference):
    chrom1, chrom2 = chrom1_chrom2_small_difference
    assert permutation_test_chrom1_lt_chrom2(chrom2, chrom1, size=100) > 0.05


def test_chrom1_gt_chrom2_no_difference(chrom1_chrom2_no_difference):
    chrom1, chrom2 = chrom1_chrom2_no_difference
    assert permutation_test_chrom1_lt_chrom2(chrom2, chrom1, size=100) > 0.05


def test_no_difference_div_by_zero(chrom1_chrom2_no_difference):
    chrom1, chrom2 = chrom1_chrom2_no_difference
    chrom2[10] = 0
    assert permutation_test_chrom1_lt_chrom2(chrom1, chrom2, size=100) > 0.05


def test_small_difference_div_by_zero(chrom1_chrom2_small_difference):
    chrom1, chrom2 = chrom1_chrom2_small_difference
    chrom2[10] = 0
    assert permutation_test_chrom1_lt_chrom2(chrom1, chrom2, size=100) <= 0.05
