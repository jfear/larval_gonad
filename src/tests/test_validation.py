import pytest

from larval_gonad.validation import GeneValidator


@pytest.mark.parametrize(
    "lit,missing,biomarkers,zscores,flag_protein,expected_score",
    [
        # Fours: Lit Gene matches BioMarkers
        ## in situ
        ({"G"}, set(), {"G"}, set(), False, 4),
        ({"G"}, set(), {"G"}, {"G"}, False, 4),
        ({"G"}, set(), {"G"}, {"G", "EPS"}, False, 4),
        ({"G", "EPS"}, set(), {"G", "EPS"}, set(), False, 4),
        ({"G", "EPS"}, set(), {"G", "EPS"}, {"C1"}, False, 4),
        ({"G", "EPS"}, {"MPS", "EC"}, {"G", "EPS", "MPS"}, set(), False, 4),
        ## protein expression
        ({"G"}, set(), {"G"}, set(), True, 4),
        ({"G", "EPS"}, set(), {"G", "EPS"}, set(), True, 4),
        # Threes: Lit Gene matches Upper Quantile of Zscores
        ## in situ
        ({"G"}, set(), set(), {"G"}, False, 3),
        ({"G", "EPS"}, set(), set(), {"G", "EPS"}, False, 3),
        ({"G", "EPS"}, {"MPS"}, set(), {"G", "EPS", "MPS"}, False, 3),
        ## protein expression
        ({"G"}, set(), set(), {"G"}, True, 3),
        ({"G", "EPS"}, set(), set(), {"G", "EPS"}, True, 3),
        # Twos: Lit Gene matches Lineage Biomarkers
        ## in situ
        ({"G", "EPS"}, set(), {"G", "EPS", "MPS"}, set(), False, 2),
        ({"C1", "C4"}, set(), {"C1", "C2", "C3"}, set(), False, 2),
        ## protein expression
        ({"EPS"}, set(), {"G"}, set(), True, 2),
        ({"C2"}, set(), {"C1"}, set(), True, 2),
        ({"C1", "C2"}, set(), {"C1", "C2", "C3"}, set(), True, 2),
        # Ones: Lit Gene matches Lineage Zscores
        ## in situ
        ({"G", "EPS"}, set(), set(), {"G", "EPS", "MPS"}, False, 1),
        ({"C1", "C4"}, set(), set(), {"C1", "C2", "C3"}, False, 1),
        ## protein expression
        ({"G", "EPS"}, set(), set(), {"G", "EPS", "MPS"}, True, 1),
        # No class
        ## in situ
        ({"G"}, set(), set(), set(), False, None),
        ## protein expression
        ({"G"}, set(), set(), set(), True, None),
        ({"MPS"}, set(), {"G"}, set(), True, None),
    ],
)
def test_GeneValidator(lit, missing, biomarkers, zscores, flag_protein, expected_score):
    gene = GeneValidator("gene1", lit, missing, biomarkers, zscores, flag_protein)
    assert gene.score == expected_score
