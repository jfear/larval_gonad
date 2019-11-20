import pytest

from larval_gonad.validation import GeneValidator


@pytest.mark.parametrize(
    "fbgn,experiment_type,lit,missing,biomarkers,zscores,expected_score",
    [
        # Fours: Lit Gene matches BioMarkers
        ("gene1", "ISH", {"G"}, set(), {"G"}, set(), 4),
        ("gene1", "IHC", {"G"}, set(), {"G"}, set(), 4),
        ("gene1", "ISH", {"G"}, set(), {"G"}, {"G"}, 4),
        ("gene1", "ISH", {"G"}, set(), {"G"}, {"G", "EPS"}, 4),
        ("gene1", "ISH", {"G", "EPS"}, set(), {"G", "EPS"}, set(), 4),
        ("gene1", "ISH", {"G", "EPS"}, set(), {"G", "EPS"}, {"C1"}, 4),
        ("gene1", "ISH", {"G", "EPS"}, {"MPS", "EC"}, {"G", "EPS", "MPS"}, set(), 4),
        # Threes: Lit Gene matches Upper Quantile of Zscores
        ("gene1", "ISH", {"G"}, set(), set(), {"G"}, 3),
        ("gene1", "IHC", {"G"}, set(), set(), {"G"}, 3),
        ("gene1", "ISH", {"G", "EPS"}, set(), set(), {"G", "EPS"}, 3),
        ("gene1", "ISH", {"G", "EPS"}, {"MPS"}, set(), {"G", "EPS", "MPS"}, 3),
        # Twos: Lit Gene matches Lineage Biomarkers
        ("gene1", "ISH", {"G", "EPS"}, set(), {"G", "EPS", "MPS"}, set(), 2),
        ("gene1", "ISH", {"C1", "C4"}, set(), {"C1", "C2", "C3"}, set(), 2),
        # Ones: Lit Gene matches Lineage Zscores
        ("gene1", "ISH", {"G", "EPS"}, set(), set(), {"G", "EPS", "MPS"}, 1),
        ("gene1", "ISH", {"C1", "C4"}, set(), set(), {"C1", "C2", "C3"}, 1),
        # No class
        ("gene1", "ISH", {"G"}, set(), set(), set(), None),
    ],
)
def test_GeneValidator(fbgn, experiment_type, lit, missing, biomarkers, zscores, expected_score):
    gene = GeneValidator(fbgn, experiment_type, lit, missing, biomarkers, zscores)
    assert gene.score == expected_score
