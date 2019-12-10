"""Summary table of ovary and testis counts.

Table has the following columns:
- FBgn
- sample_ID
- Count
- stage = 'L3'
- tissue
- rep
- chrom

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
    )


def tissue_mapper(sample_ID: pd.Series):
    mapper = dict(TCP="testis", OCP="ovary")
    return sample_ID.str.extract(r".*_(\w+)", expand=False).map(mapper)


def rep_mapper(sample_ID: pd.Series):
    mapper = {
        "A1": "1",
        "A2": "2",
        "A3": "3",
        "A4": "4",
        "A5": "1",
        "A6": "2",
        "A7": "3",
        "A8": "4",
    }
    return sample_ID.str.extract(r"(.*)_\w+", expand=False).map(mapper)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="bulk2-rnaseq-wf",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                counts="../output/bulk2-rnaseq-wf/rnaseq_aggregation/gene_level_counts.tsv",
            ),
        )

    main()
