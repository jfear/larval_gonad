"""Generate a table mapping chromosome locations for different species.

1. This script uses the ortholog tale to map YOgn IDs to FBgn IDs.
2. It then uses the gtf to figure out chromosome locations for each gene
   across species.

"""
import re

import pandas as pd

from larval_gonad.io import GffRow


def main():
    # create dict mapper species to ortholog tsv
    ortho_mapper = dict(
        ((re.findall(r"(d\w+)\.tsv", fname)[0], fname) for fname in snakemake.input.orthologs)
    )

    # create dict mapper species to gtf
    gtf_mapper = dict(
        ((re.findall(r"(d\w+)\.gtf", fname)[0], fname) for fname in snakemake.input.gtf)
    )

    # get chromosome locations for each gene
    locations = pd.concat(
        (
            get_chromosome_locations(species, gtf_mapper, ortho_mapper)
            for species in gtf_mapper.keys()
        ),
        sort=True,
        axis=1,
    )

    locations.rename_axis("FBgn").reset_index().to_feather(snakemake.output[0])


def get_chromosome_locations(species, gtf_mapper, ortho_mapper):
    """Chromosome locations for each gene for the provided species.

    1. Get all gene locations from GTF.
    2. Use the ortholog table to convert YOgn ids to FBgn

    """
    gtf_fname = gtf_mapper[species]
    ortho_fname = ortho_mapper[species]

    # Get gene chrom locations for species
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
    DEBUG = False

    if DEBUG:
        from larval_gonad.debug import snakemake_debug
        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input={
                "orthologs": [
                    "../output/expression-atlas-wf/orthologs/dana.tsv",
                    "../output/expression-atlas-wf/orthologs/dmel.tsv",
                    "../output/expression-atlas-wf/orthologs/dmoj.tsv",
                    "../output/expression-atlas-wf/orthologs/dper.tsv",
                    "../output/expression-atlas-wf/orthologs/dpse.tsv",
                    "../output/expression-atlas-wf/orthologs/dvir.tsv",
                    "../output/expression-atlas-wf/orthologs/dwil.tsv",
                    "../output/expression-atlas-wf/orthologs/dyak.tsv",
                ],
                "gtf": [
                    "../output/expression-atlas-wf/GTF/dana.gtf",
                    "../output/expression-atlas-wf/GTF/dmel.gtf",
                    "../output/expression-atlas-wf/GTF/dmoj.gtf",
                    "../output/expression-atlas-wf/GTF/dper.gtf",
                    "../output/expression-atlas-wf/GTF/dpse.gtf",
                    "../output/expression-atlas-wf/GTF/dvir.gtf",
                    "../output/expression-atlas-wf/GTF/dwil.gtf",
                    "../output/expression-atlas-wf/GTF/dyak.gtf",
                ],
            }

        )
    main()
