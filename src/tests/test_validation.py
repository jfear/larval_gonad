import pytest

from larval_gonad.validation import GeneValidator


@pytest.mark.parametrize(
    "fbgn,lit,biomarkers,zscores,expected_score",
    [
        ("gene1", {"G"}, {"G"}, {}, 4),
        ("gene1", {"G"}, {"G"}, {"G"}, 4),
        ("gene1", {"G"}, {"G"}, {"G", "EPS"}, 4),
        ("gene1", {"G", "EPS"}, {"G", "EPS"}, {}, 4),
        ("gene1", {"G", "EPS"}, {"G", "EPS"}, {"C1"}, 4),
    ],
)
def test_GeneValidator(fbgn, lit, biomarkers, zscores, expected_score):
    gene = GeneValidator(fbgn, lit, biomarkers, zscores)
    assert gene.score == expected_score
