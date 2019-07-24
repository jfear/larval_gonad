import re

from lcdblib.utils.chrom_convert import import_conversion
from larval_gonad.io import GffRow

GTF = snakemake.input[0]
OUTNAME = snakemake.output[0]


def main():
    fb_mapper = import_conversion("UCSC", "FlyBase")

    with open(GTF) as fh, open(OUTNAME, "w") as fo:
        fo.write(
            "FBgn\tgene_symbol\tUCSC_chrom\tFB_chrom\tstart\tend\tlength\tstrand\n"
        )
        for row in fh.readlines():
            grow = GffRow(row)
            if grow.is_gene:
                fo.write(
                    "\t".join(
                        [
                            grow["gene_id"],
                            grow["gene_symbol"],
                            grow.seqid,
                            fb_mapper[grow.seqid],
                            grow.start,
                            grow.end,
                            str(int(grow.end) - int(grow.start)),
                            grow.strand,
                        ]
                    )
                    + "\n"
                )
