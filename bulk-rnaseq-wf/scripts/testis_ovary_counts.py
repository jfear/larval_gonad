"""Summary table of ovary and testis counts.

Table has the following columns:
- FBgn
- sample_ID
- Count
- stage = 'L3'
- tissue
- rep
- chrom
- data_source = 'RNA-Seq'

"""
import os
import pandas as pd


def main():
    fbgn2chrom = get_fbgn2chrom()
    counts = get_counts()
    df = counts.join(fbgn2chrom, how="inner").reset_index()
    df.to_feather(snakemake.output[0])


def get_fbgn2chrom():
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .query('FB_chrom in ["X", "2L", "2R", "3L", "3R", "4", "Y"]')
        .squeeze()
        .rename("chrom")
    )


def get_counts():
    return (
        pd.read_csv(snakemake.input.counts, sep="\t")
        .pipe(lambda x: x[~x.Geneid.str.contains("ERCC")])
        .melt(id_vars="Geneid", var_name="sample_ID", value_name="Count")
        .set_index("Geneid")
        .rename_axis("FBgn")
        .assign(stage="L3")
        .assign(tissue=lambda x: tissue_mapper(x.sample_ID))
        .assign(rep=lambda x: rep_mapper(x.sample_ID))
        .assign(data_source="RNA-Seq")
    )


def tissue_mapper(sample_ID: pd.Series):
    mapper = dict(TCP="testis", OCP="ovary")
    return sample_ID.str.extract(r".*_(\w+)", expand=False).map(mapper)


def rep_mapper(sample_ID: pd.Series):
    mapper = {
        # Testis
        "B5": "1",
        "B6": "2",
        "B7": "3",
        "B8": "4",
        "A5": "5",
        "A6": "6",
        "A7": "7",
        "A8": "8",
        # Ovary
        "B9": "1",
        "B10": "2",
        "B11": "3",
        "B12": "4",
        "A1": "5",
        "A2": "6",
        "A3": "7",
        "A4": "8",
    }
    return sample_ID.str.extract(r"(.*)_\w+", expand=False).map(mapper)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="bulk-rnaseq-wf",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                counts="../output/bulk-rnaseq-wf/rnaseq_aggregation/gene_level_counts.tsv",
            ),
        )

    main()
