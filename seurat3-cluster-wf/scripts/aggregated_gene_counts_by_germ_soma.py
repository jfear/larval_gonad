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
        pd.read_feather(snakemake.input.counts)
        .assign(cell_type=lambda x: cell_class_mapper(x.cluster))
        .groupby(["FBgn", "rep", "cell_type"])
        .UMI.sum()
        .rename("Count")
        .reset_index()
        .set_index("FBgn")
        .assign(tissue="testis")
        .assign(stage="L3")
        .assign(sample_ID=lambda x: "L3_scRNAseq_" + x.rep)
        .assign(data_source="scRNA-Seq")
    )

def cell_class_mapper(cluster):
    mapper = {
        "G": "Germline",
        "EPS": "Germline",
        "MPS": "Germline",
        "LPS": "Germline",
        "C1": "Cyst Lineage",
        "C2": "Cyst Lineage",
        "C3": "Cyst Lineage",
        "C4": "Cyst Lineage",
        "T": "Other Somatic",
        "P": "Other Somatic",
    }
    return cluster.map(mapper)

if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                counts="../output/seurat3-cluster-wf/raw_by_cluster_rep.feather",
            ),
        )

    main()
