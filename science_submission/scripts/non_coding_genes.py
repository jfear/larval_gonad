from pickle import dump

from gffutils import FeatureDB

fname = snakemake.input[0]
oname = snakemake.output[0]


def main():
    db = FeatureDB(fname)

    noncoding = set()
    for gene in db.features_of_type('gene'):
        if len(list(db.children(gene, featuretype='CDS'))) == 0:
            noncoding.add(gene.id)

    with open(oname, 'wb') as fh:
        dump(list(noncoding), fh)


if __name__ == '__main__':
    main()
