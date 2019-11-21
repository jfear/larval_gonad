"""Cell type validation logic. """
from pathlib import Path
from itertools import product

import numpy as np

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
        if self.lit_gene == set():
            return

        if self.flag_protein:
            if self.lvl1_protein():
                return
            elif self.lvl2_protein():
                return
            elif self.lvl3_protein():
                return
            elif self.lvl4_protein():
                return
        else:
            if self.lvl1_rna():
                return
            elif self.lvl2_rna():
                return
            elif self.lvl3_rna():
                return
            elif self.lvl4_rna():
                return

    def lvl1_rna(self):
        """{Lit} == {Biomarkers}"""
        if self.lit_gene == self.non_missing_biomarker:
            self.score = 4
            return True

    def lvl1_protein(self):
        """{Lit} == {Biomarkers} or {Biomarkers} less than"""
        if self.lit_gene == self.non_missing_biomarker:
            self.score = 4
            return True
        elif self.all_proceeds(self.lit_gene, self.biomarker, self.germ_lineage):
            self.score = 4
            return True
        elif self.all_proceeds(self.lit_gene, self.biomarker, self.cyst_lineage):
            self.score = 4
            return True

    def lvl2_rna(self):
        """{Lit} == {q75(Zscores)}"""
        if self.lit_gene == self.non_missing_zscore:
            self.score = 3
            return True

    def lvl2_protein(self):
        """{Lit} == {q75(Zscores)} or {q75(Zscores)} less than"""
        if self.lit_gene == self.non_missing_zscore:
            self.score = 3
            return True
        elif self.all_proceeds(self.lit_gene, self.zscore, self.germ_lineage):
            self.score = 3
            return True
        elif self.all_proceeds(self.lit_gene, self.zscore, self.cyst_lineage):
            self.score = 3
            return True

    def lvl3_rna(self):
        """({Lit} subset {Germ} & {Biomarkers} subset {Germ}) or
           ({Lit} subset {Cyst} & {Biomarkers} subset {Cyst})
        """
        # Germ line
        lit_lineage = self.lit_gene.intersection(self.germ_lineage)
        biomarker_lineage = self.biomarker.intersection(self.germ_lineage)

        if (lit_lineage != set()) & (biomarker_lineage != set()):
            self.score = 2
            return True

        # Cyst line
        lit_lineage = self.lit_gene.intersection(self.cyst_lineage)
        biomarker_lineage = self.biomarker.intersection(self.cyst_lineage)

        if (lit_lineage != set()) & (biomarker_lineage != set()):
            self.score = 2
            return True

    def lvl3_protein(self):
        """Gene expression proceeds in the lineage."""
        if self.any_proceeds(self.lit_gene, self.biomarker, self.germ_lineage):
            self.score = 2
            return True

        # Cyst line
        lit_lineage = self.lit_gene.intersection(self.cyst_lineage)
        biomarker_lineage = self.biomarker.intersection(self.cyst_lineage)

        if (lit_lineage != set()) & (biomarker_lineage != set()):
            self.score = 2
            return True

    def lvl4_rna(self):
        """({Lit} subset {Germ} & {q75(Zscores)} subset {Germ}) or
           ({Lit} subset {Cyst} & {q75(Zscores)} subset {Cyst})
        """
        # Germ line
        lit_lineage = self.lit_gene.intersection(self.germ_lineage)
        zscore_lineage = self.zscore.intersection(self.germ_lineage)

        if (lit_lineage != set()) & (zscore_lineage != set()):
            self.score = 1
            return True

        # Cyst line
        lit_lineage = self.lit_gene.intersection(self.cyst_lineage)
        zscore_lineage = self.zscore.intersection(self.cyst_lineage)

        if (lit_lineage != set()) & (zscore_lineage != set()):
            self.score = 1
            return True

    def lvl4_protein(self):
        """Gene expression any_proceeds in the lineage."""
        if self.any_proceeds(self.lit_gene, self.zscore, self.germ_lineage):
            self.score = 1
            return True

        # Cyst line
        lit_lineage = self.lit_gene.intersection(self.cyst_lineage)
        zscore_lineage = self.zscore.intersection(self.cyst_lineage)

        if (lit_lineage != set()) & (zscore_lineage != set()):
            self.score = 1
            return True

    def rna_proceeds_any(self, rna, proteins, lineage):
        for protein in proteins:
            protein_idx = lineage.index(protein) + 1
            rna_idx = lineage.index(rna) + 1
            if rna_idx <= protein_idx:
                return True
        return False

    def protein_proceeds_any(self, rnas, protein, lineage):
        for rna in rnas:
            protein_idx = lineage.index(protein) + 1
            rna_idx = lineage.index(rna) + 1
            if protein_idx <= rna_idx:
                return True
        return False

    def all_proceeds(self, protein_cell_types, cell_types, lineage):
        """Check if all cell types proceed or are the same as the protein index."""
        protein_lineage = protein_cell_types.intersection(lineage)
        rna_lineage = cell_types.intersection(lineage)

        if (protein_lineage == set()) | (rna_lineage == set()):
            return 

        _rna_before_protein = np.array([
            self.rna_proceeds_any(rna, protein_lineage, lineage)
            for rna in rna_lineage
        ])

        _protein_before_rna = np.array([
            self.protein_proceeds_any(rna_lineage, protein, lineage)
            for protein in protein_lineage
        ])

        if np.all(_rna_before_protein) and not np.all(_protein_before_rna):
            return True

    def any_proceeds(self, protein_cell_types, cell_types, lineage):
        """Check if any cell types proceed the protein index."""
        protein_lineage = protein_cell_types.intersection(lineage)
        rna_lineage = cell_types.intersection(lineage)

        if (protein_lineage == set()) | (rna_lineage == set()):
            return 

        _rna_before_protein = np.array([
            self.rna_proceeds_any(rna, protein_lineage, lineage)
            for rna in rna_lineage
        ])

        _protein_before_rna = np.array([
            self.protein_proceeds_any(rna_lineage, protein, lineage)
            for protein in protein_lineage
        ])

        if np.any(_rna_before_protein):
            return True

    def __str__(self):
        return f"{self.fbgn}\t{self.lit_gene}\t{self.biomarker}\t{self.zscore}\t{self.score}"
