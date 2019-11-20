"""Cell type validation logic. """
from pathlib import Path

from larval_gonad.config import read_config


config = read_config(Path(Path(__file__).parent, "../../config/common.yaml"))


class GeneValidator(object):
    germ_lineage = config["germ"]
    cyst_lineage = config["cyst"]

    def __init__(self, FBgn, lit_gene, lit_missing, biomarker, zscore, flag_protein=False):
        self.fbgn = FBgn
        self.lit_gene = lit_gene
        self.lit_missing = lit_missing
        self.biomarker = biomarker
        self.non_missing_biomarker = self.biomarker - self.lit_missing
        self.zscore = zscore
        self.flag_protein = flag_protein
        self.non_missing_zscore = self.zscore - self.lit_missing
        self.score = None
        self.validate()

    def validate(self):
        if self.lvl1():
            return
        elif self.lvl2():
            return

        if self.flag_protein:
            if self.lvl3_protein():
                return
            elif self.lvl4_protein():
                return
        else:
            if self.lvl3_rna():
                return
            elif self.lvl4_rna():
                return

    def lvl1(self):
        """{Lit} == {Biomarkers}"""
        if self.lit_gene == self.non_missing_biomarker:
            self.score = 4
            return True

    def lvl2(self):
        """{Lit} == {q75(Zscores)}"""
        if self.lit_gene == self.non_missing_zscore:
            self.score = 3
            return True

    def lvl3_rna(self):
        """({Lit} subset {Germ} & {Biomarkers} subset {Germ}) or
           ({Lit} subset {Cyst} & {Biomarkers} subset {Cyst})
        """
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

    def lvl3_protein(self):
        """Gene expression proceeds in the lineage."""
        if (
            self.lit_gene.issubset(self.germ_lineage)
            & self.biomarker.issubset(self.germ_lineage)
            & bool(self.biomarker)
        ):
            protein_idx = max([self.germ_lineage.index(x) + 1 for x in self.lit_gene])
            if self.proceeds(protein_idx, self.biomarker, self.germ_lineage):
                self.score = 2
                return True
        elif (
            self.lit_gene.issubset(self.cyst_lineage)
            & self.biomarker.issubset(self.cyst_lineage)
            & bool(self.biomarker)
        ):
            protein_idx = max([self.cyst_lineage.index(x) + 1 for x in self.lit_gene])
            if self.proceeds(protein_idx, self.biomarker, self.cyst_lineage):
                self.score = 2
                return True

    def lvl4_rna(self):
        """({Lit} subset {Germ} & {q75(Zscores)} subset {Germ}) or
           ({Lit} subset {Cyst} & {q75(Zscores)} subset {Cyst})
        """
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

    def lvl4_protein(self):
        """Gene expression proceeds in the lineage."""
        if (
            self.lit_gene.issubset(self.germ_lineage)
            & self.zscore.issubset(self.germ_lineage)
            & bool(self.zscore)
        ):
            protein_idx = max([self.germ_lineage.index(x) + 1 for x in self.lit_gene])
            if self.proceeds(protein_idx, self.zscore, self.germ_lineage):
                self.score = 1
                return True
        elif (
            self.lit_gene.issubset(self.cyst_lineage)
            & self.zscore.issubset(self.cyst_lineage)
            & bool(self.zscore)
        ):
            protein_idx = max([self.cyst_lineage.index(x) + 1 for x in self.lit_gene])
            if self.proceeds(protein_idx, self.zscore, self.cyst_lineage):
                self.score = 1
                return True

    def proceeds(self, protein_idx, cell_types, lineage):
        """Check if any cell types proceed the protein index."""
        for cell_type in cell_types:
            rna_idx = lineage.index(cell_type) + 1
            _diff = protein_idx - rna_idx
            if _diff in [0, 1]:
                return True

    def __str__(self):
        return f"{self.fbgn}\t{self.lit_gene}\t{self.biomarker}\t{self.zscore}\t{self.score}"
