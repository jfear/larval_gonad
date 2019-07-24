from pickle import dump
import pandas as pd

INPUT_FILE = snakemake.params[0]
OUTPUT_FILE = snakemake.output[0]


def main():
    df = pd.read_csv(
        INPUT_FILE,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "gene_symbol",
            "organism",
            "primary_FBgn",
            "secondary_FBgn",
            "annotation_ID",
            "secondary_annotation_ID",
        ],
    ).set_index("primary_FBgn")

    mapper = {}
    for idx, row in df.iterrows():
        mapper[idx] = idx
        try:
            for secondary in row.secondary_FBgn.split(","):
                mapper[secondary] = idx
        except AttributeError:
            pass

    with open(OUTPUT_FILE, 'wb') as fh:
        dump(mapper, fh)


if __name__ == "__main__":
    main()
