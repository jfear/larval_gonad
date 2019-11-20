
"""Cell type validation logic. """

class GeneValidator(object):
    def __init__(self, FBgn, lit_gene, biomarkers, zscores):
        self.fbgn = FBgn
        self.lit_gene = lit_gene
        self.biomarker = biomarkers.get(self.fbgn, set())
        self.zscore = zscores.get(self.fbgn, set())
        self.score = None
        self.validate()

    def validate(self):
        if self.lvl1():
            return

        if self.lvl2():
            return

        if self.lvl3():
            return

    def lvl1(self):
        """Lit and Biomarkers exact match"""
        if self.lit_gene == self.biomarker:
            self.score = 4
            return True

    def lvl2(self):
        """Lit and upper zscore quantail exact match"""
        if self.lit_gene == self.zscore:
            self.score = 3
            return True

    def lvl3(self):
        pass

    def __str__(self):
        return f"{self.fbgn}\t{self.lit_gene}\t{self.biomarker}\t{self.zscore}\t{self.score}"
