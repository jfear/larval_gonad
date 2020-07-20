import pytest
import numpy as np
import pandas as pd

from larval_gonad.stats import PairwisePermutationTest, permute_sample

np.random.seed(42)

SAMPLE_SIZE = 1_000
DIFFERENCES = {"large": 500, "small": 50, "minute": 5, "none": 0}


@pytest.fixture(scope="session")
def toy_data():
    data1 = np.array([1, 2, 3, 4, 5])
    data2 = np.array([11, 12, 13, 14])
    return data1, data2


@pytest.fixture(scope="session")
def random_counts():
    return np.random.lognormal(6, size=SAMPLE_SIZE)


@pytest.fixture(scope="session", params=DIFFERENCES.values(), ids=DIFFERENCES.keys())
def data1_data2(random_counts, request):
    data1 = random_counts
    data2 = data1 + np.random.normal(request.param, size=SAMPLE_SIZE)
    return data1, data2


@pytest.fixture(scope="session")
def dataframe(data1_data2):
    data1, data2 = data1_data2
    data = np.concatenate((data1, data2))
    labels = ["data1"] * len(data1) + ["data2"] * len(data2)
    return pd.DataFrame({"data": data, "labels": labels})


def test_permute_sample(toy_data):
    data1, data2 = toy_data
    permuted_sample_1, permuted_sample_2 = permute_sample(
        data1, data2, seed=42
    )
    assert len(permuted_sample_1) == len(data1)
    assert np.all(permuted_sample_1 == np.array([13, 2, 11, 1, 14]))
    assert len(permuted_sample_2) == len(data2)
    assert np.all(permuted_sample_2 == np.array([3, 5, 4, 12]))


def test_differences(dataframe):
    pw = PairwisePermutationTest(
        group_column="labels", value_column="data", data=dataframe, n_permutations=10
    ).fit()
    res = pw._results[0]
    if np.abs(res.log2FoldChange) <= 0.02:
        # Really small difference, expect non-significant p-value
        assert res.p_value >= 0.05
    else:
        # Bigger difference, expect significant p-value
        assert res.p_value <= 0.05
