import pytest

from larval_gonad.validation import GeneValidator


@pytest.mark.parametrize(
    "lit,missing,biomarkers,zscores,flag_protein,expected_score",
    [
        # Fours: Lit Gene matches BioMarkers
        ## in situ
        ({"G"}, set(), {"G"}, set(), False, 4),
        ({"G", "EPS"}, set(), {"G", "EPS"}, set(), False, 4),
        ({"G", "EPS"}, {"MPS", "EC"}, {"G", "EPS", "MPS"}, set(), False, 4),
        ## protein expression
        ({"G"}, set(), {"G"}, set(), True, 4),
        ({"EPS"}, set(), {"G"}, set(), True, 4),
        ({"G", "EPS"}, set(), {"G", "EPS"}, set(), True, 4),
        ({"EPS", "MPS", "LPS"}, set(), {"G", "EPS"}, set(), True, 4),
        ({"C2"}, set(), {"C1"}, set(), True, 4),
        ({"C1", "C2", "C3"}, set(), {"C1", "C2"}, set(), True, 4),

        # Threes: Lit Gene matches Upper Quantile of Zscores
        ## in situ
        ({"G"}, set(), set(), {"G"}, False, 3),
        ({"G", "EPS"}, set(), set(), {"G", "EPS"}, False, 3),
        ({"G", "EPS"}, {"MPS", "EC"}, set(), {"G", "EPS", "MPS"}, False, 3),
        ## protein expression
        ({"G"}, set(), set(), {"G"}, True, 3),
        ({"EPS"}, set(), set(), {"G"}, True, 3),
        ({"G", "EPS"}, set(), set(), {"G", "EPS"}, True, 3),
        ({"EPS", "MPS", "LPS"}, set(), set(), {"G", "EPS"}, True, 3),
        ({"C2"}, set(), set(), {"C1"}, True, 3),
        ({"C1", "C2", "C3"}, set(), set(), {"C1", "C2"}, True, 3),

        # Twos: Lit Gene matches Lineage Biomarkers
        ## in situ
        ({"G", "EPS"}, set(), {"G", "EPS", "MPS"}, set(), False, 2),
        ({"C1", "C4"}, set(), {"C1", "C2", "C3"}, set(), False, 2),
        ({"EPS", "MPS", "LPS", "P"}, set(), {"EPS"},set(),  False, 2),

        ## protein expression
        ({"C1", "C2"}, set(), {"C1", "C2", "C3"}, set(), True, 2),
        ({"G", "EPS"}, set(), {"EPS", "MPS"}, set(), True, 2),
        ({"G", "EPS", "MPS", "P"}, set(), {"EPS", "MPS"}, set(), True, 2),
        ({"G", "EPS", "MPS"}, set(), {"EPS", "MPS", "P"}, set(), True, 2),

        # Ones: Lit Gene matches Lineage Zscores
        ## in situ
        ({"G", "EPS"}, set(), set(), {"G", "EPS", "MPS"}, False, 1),
        ({"C1", "C4"}, set(), set(), {"C1", "C2", "C3"}, False, 1),
        ({"EPS", "MPS", "LPS", "P"}, set(), set(), {"EPS"}, False, 1),
        ## protein expression
        ({"G", "EPS"}, set(), set(), {"G", "EPS", "MPS"}, True, 1),

        # No class
        (set(), {"G"}, set(), set(), False, None),
        ## in situ
        ({"G"}, set(), set(), set(), False, None),
        ## protein expression
        ({"G"}, set(), set(), set(), True, None),
    ],
)
def test_GeneValidator(lit, missing, biomarkers, zscores, flag_protein, expected_score):
    gene = GeneValidator("gene1", lit, missing, biomarkers, zscores, flag_protein)
    assert gene.score == expected_score
