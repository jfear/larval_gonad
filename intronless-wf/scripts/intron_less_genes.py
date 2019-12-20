from pickle import dump

from gffutils import FeatureDB

fname = snakemake.input[0]
oname = snakemake.output[0]


def main():
    db = FeatureDB(fname)

    intron_counts = {}
    for gene in db.features_of_type('gene'):
        intron_counts[gene.id] = 0

    for intron in db.create_introns():
        intron_counts[intron.attributes['gene_id'][0]] += 1

    intron_less = []
    for fbgn, cnt in intron_counts.items():
        if cnt == 0:
            intron_less.append(fbgn)

    with open(oname, 'wb') as fh:
        dump(intron_less, fh)


if __name__ == '__main__':
    main()
