import re

import pandas as pd

from larval_gonad.io import GffRow

ORTHROLOG_FILES = snakemake.input.orthologs
GTF_FILES = snakemake.input.gtf
OUTPUT_FILE = snakemake.output[0]

# Debug settings
# import os
# os.chdir('neox-wf/scripts')
# ORTHROLOG_FILES = [
#     '../../output/neox-wf/orthologs/dana.tsv',
#     '../../output/neox-wf/orthologs/dmel.tsv',
#     '../../output/neox-wf/orthologs/dmoj.tsv',
#     '../../output/neox-wf/orthologs/dper.tsv',
#     '../../output/neox-wf/orthologs/dpse.tsv',
#     '../../output/neox-wf/orthologs/dvir.tsv',
#     '../../output/neox-wf/orthologs/dwil.tsv',
#     '../../output/neox-wf/orthologs/dyak.tsv'
# ]
# GTF_FILES = [
#     '../../output/neox-wf/GTF/dana.gtf',
#     '../../output/neox-wf/GTF/dmel.gtf',
#     '../../output/neox-wf/GTF/dmoj.gtf',
#     '../../output/neox-wf/GTF/dper.gtf',
#     '../../output/neox-wf/GTF/dpse.gtf',
#     '../../output/neox-wf/GTF/dvir.gtf',
#     '../../output/neox-wf/GTF/dwil.gtf',
#     '../../output/neox-wf/GTF/dyak.gtf'
# ]


def main():
    ortho_mapper = dict(
        ((re.findall(r"(d\w+)\.tsv", fname)[0], fname) for fname in ORTHROLOG_FILES)
    )

    gtf_mapper = dict(((re.findall(r"(d\w+)\.gtf", fname)[0], fname) for fname in GTF_FILES))

    locations = pd.concat((
        get_locations(species, gtf_mapper, ortho_mapper)
        for species in gtf_mapper.keys()
    ), sort=True, axis=1)

    locations.rename_axis('FBgn').reset_index().to_feather(OUTPUT_FILE)


def get_locations(species, gtf_mapper, ortho_mapper):
    gtf_fname = gtf_mapper[species]
    ortho_fname = ortho_mapper[species]

    # Get gene locations for species
    results = set()
    with open(gtf_fname) as fh:
        for row in fh.read().strip().split("\n"):
            gtf = GffRow(row)
            results.add((gtf.parsed_attributes["gene_id"], gtf.seqid))
    locations = pd.DataFrame(results, columns=["YOgnID", "chrom"]).set_index("YOgnID")

    # Convert YOgnIDS to FBgn
    yo2FBgn = (
        pd.read_csv(ortho_fname, sep="\t", usecols=["YOgnID", "Dmel"], index_col="YOgnID")
        .squeeze()
        .rename("FBgn")
    )

    return locations.join(yo2FBgn).dropna().set_index("FBgn", drop=True).squeeze().rename(species)


if __name__ == "__main__":
    main()
