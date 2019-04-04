import os
from pathlib import Path

import pandas as pd
import gffutils


def main():
    oname = Path('../output/gene_ts_lengths.tsv')
    gtf = Path(os.environ['REFERENCES_DIR'], 'dmel/r6-16/gtf/dmel_r6-16.gtf.db')

    db = gffutils.FeatureDB(gtf.as_posix())

    gene_ts_lengths = []
    for gene in db.features_of_type('gene'):
        length = 0
        for exon in db.merge(db.children(gene, featuretype='exon'), ignore_strand=True):
            length += len(exon)

        try:
            assert len(gene) >= length
            gene_ts_lengths.append([gene.id, length])
        except AssertionError:
            print(gene.id, len(gene), length)

    gene_ts_lengths = pd.DataFrame(gene_ts_lengths, columns=['FBgn', 'gene_ts_length']).set_index('FBgn')
    gene_ts_lengths.to_csv(oname, sep='\t')


if __name__ == '__main__':
    main()
