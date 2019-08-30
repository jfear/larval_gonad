"""Parse GTF file and aggregate gene metadata.

Capture key metadata from GTF file:
    * chromosome/scaffold
    * start
    * end
    * strand
"""
import pandas as pd

from larval_gonad.io import GffRow


def main():
    res = [
        (ts.parsed_attributes["gene_id"], ts.seqid, ts.start, ts.end, ts.strand)
        for ts in gff_reader(snakemake.input[0])
    ]

    df = pd.DataFrame(res, columns=["YOgn", "chrom", "start", "end", "strand"])

    yo2metdata = df.groupby("YOgn").agg(
        {"chrom": "first", "start": min, "end": max, "strand": "first"}
    )

    yo2metdata.reset_index().to_feather(snakemake.output[0])


def gff_reader(file_name):
    with open(file_name) as fh:
        for row in fh:
            gff = GffRow(row)
            if gff.type == "transcript":
                yield gff


if __name__ == "__main__":
    DEBUG = False

    if DEBUG:
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf", input="../output/expression-atlas-wf/GTF/dana.gtf"
        )

    main()
