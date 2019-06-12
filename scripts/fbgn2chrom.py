from gff import GffRow

GTF = snakemake.input[0]
OUTNAME = snakemake.output[0]


def main():
    with open(GTF) as fh, open(OUTNAME, "w") as fo:
        fo.write("FBgn\tchrom\n")
        for row in fh.readlines():
            grow = GffRow(row)
            if grow.is_gene:
                fo.write("\t".join([grow["gene_id"], grow.seqid]) + "\n")


if __name__ == "__main__":
    main()
