import pandas as pd

URL = snakemake.params[0]
OUTPUT_FILE = snakemake.output[0]
COLUMN_ORDER = ["FBgn", "dana", "dmoj", "dper", "dpse", "dvir", "dwil", "dyak", "dmel"]


def main():
    col_names = [
        "FBgn",
        "gene_symbol",
        "chrom",
        "location",
        "strand",
        "ortholog_FBgn",
        "ortholog_gene_symbol",
        "ortholog_chrom",
        "ortholog_location",
        "ortholog_strand",
        "ortholog_group_id",
    ]
    ortholog_table = pd.read_csv(URL, header=None, sep="\t", names=col_names, comment="#").assign(
        species=lambda x: x.ortholog_gene_symbol.str.extract(r"(\w+)\\.*", expand=False).str.lower()
    )

    dmel_chroms = ortholog_table.set_index("FBgn").chrom.rename("dmel")

    pivoted = pd.pivot_table(
        ortholog_table, index="FBgn", columns="species", values="ortholog_chrom", aggfunc="first"
    )

    fb_orthologs = (
        pivoted.join(dmel_chroms).reset_index()[COLUMN_ORDER].drop_duplicates().reset_index()
    )

    fb_orthologs.to_feather(OUTPUT_FILE)


if __name__ == "__main__":
    main()
