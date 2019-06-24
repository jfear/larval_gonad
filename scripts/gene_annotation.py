import re

from lcdblib.utils.chrom_convert import import_conversion

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


class GffRow(object):
    def __init__(self, row):
        self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes = row.strip().split(
            "\t"
        )
        self.is_gene = self.type == "gene"
        self.parsed_attributes = self.parse_attributes()

    def parse_attributes(self):
        parsed_attributes = {}
        for attr in self.attributes.split(";"):
            mm = re.search('(?P<key>.*?)\s+"(?P<value>.*?)"', attr)
            if mm:
                parsed_attributes[mm.group("key").strip()] = mm.group("value").strip()
        return parsed_attributes

    def __getitem__(self, key):
        return self.parsed_attributes[key]


if __name__ == "__main__":
    main()
