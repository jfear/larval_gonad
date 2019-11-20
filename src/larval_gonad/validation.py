"""Cell type validation logic. """
from pathlib import Path

from larval_gonad.config import read_config


config = read_config(Path(Path(__file__).parent, "../../config/common.yaml"))


class GeneValidator(object):
    germ_lineage = config["germ"]
    cyst_lineage = config["cyst"]

    def __init__(self, FBgn, experiment_type, lit_gene, lit_missing, biomarker, zscore):
        self.fbgn = FBgn
        self.exp_type = experiment_type
        self.lit_gene = lit_gene
        self.lit_missing = lit_missing
        self.biomarker = biomarker
        self.non_missing_biomarker = self.biomarker - self.lit_missing
        self.zscore = zscore
        self.non_missing_zscore = self.zscore - self.lit_missing
        self.score = None
        self.validate()

    def validate(self):
        if self.lvl1():
            return
        elif self.lvl2():
            return

        if self.exp_type == "ISH":
            if self.lvl3_ish():
                return
            elif self.lvl4_ish():
                return

        if self.exp_type == "IHC":
            if self.lvl3_ihc():
                return
            elif self.lvl4_ihc():
                return

    def lvl1(self):
        """Lit and Biomarkers exact match"""
        if self.lit_gene == self.non_missing_biomarker:
            self.score = 4
            return True

    def lvl2(self):
        """Lit and upper zscore quantail exact match"""
        if self.lit_gene == self.non_missing_zscore:
            self.score = 3
            return True

    def lvl3_ish(self):
        if (
            self.lit_gene.issubset(self.germ_lineage)
            & self.biomarker.issubset(self.germ_lineage)
            & bool(self.biomarker)
        ) | (
            self.lit_gene.issubset(self.cyst_lineage)
            & self.biomarker.issubset(self.cyst_lineage)
            & bool(self.biomarker)
        ):
            self.score = 2
            return True

    def lvl3_ihc(self):
        if (
            self.lit_gene.issubset(self.germ_lineage)
            & self.biomarker.issubset(self.germ_lineage)
            & bool(self.biomarker)
        ) | (
            self.lit_gene.issubset(self.cyst_lineage)
            & self.biomarker.issubset(self.cyst_lineage)
            & bool(self.biomarker)
        ):
            self.score = 2
            return True

    def lvl4_ish(self):
        if (
            self.lit_gene.issubset(self.germ_lineage)
            & self.zscore.issubset(self.germ_lineage)
            & bool(self.zscore)
        ) | (
            self.lit_gene.issubset(self.cyst_lineage)
            & self.zscore.issubset(self.cyst_lineage)
            & bool(self.zscore)
        ):
            self.score = 1
            return True

    def lvl4_ihc(self):
        if (
            self.lit_gene.issubset(self.germ_lineage)
            & self.zscore.issubset(self.germ_lineage)
            & bool(self.zscore)
        ) | (
            self.lit_gene.issubset(self.cyst_lineage)
            & self.zscore.issubset(self.cyst_lineage)
            & bool(self.zscore)
        ):
            self.score = 1
            return True

    def __str__(self):
        return f"{self.fbgn}\t{self.lit_gene}\t{self.biomarker}\t{self.zscore}\t{self.score}"
