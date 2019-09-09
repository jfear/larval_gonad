library(Seurat)
library(dplyr)
library(tibble)
library(feather)
library(ggplot2)

if (Sys.getenv("SNAKE_DEBUG", FALSE)) {
    source("../lib/snakemake.R")
    snakemake <- snakemake_debug(
        input = list(
            robj = "../output/seurat3-cluster-wf/combined_n3.Robj",
            gene_annot = "../references/gene_annotation_dmel_r6-26.feather"
        ),
        params = list(germ_clusters = c(6, 4, 2, 0))
    )
}

# Load Metadata
fbgn2symbol <- read_feather(snakemake@input[["gene_annot"]], columns = c("FBgn", "gene_symbol"))

# Load combined robj
load(snakemake@input[["robj"]])

# Subset out the germline clusters and clear metadata
sobj <- subset(combined, ident = snakemake@params[["germ_clusters"]])
sobj@meta.data <- sobj@meta.data %>%
    select(orig.ident, nCount_RNA, nFeature_RNA, sample)

# Re-process the subset
sobj <- RunPCA(sobj, npcs = 30, features = VariableFeatures(sobj), verbose = FALSE)
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:9)
sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:9, verbose = FALSE)
sobj <- FindClusters(sobj, resolution = 0.4, verbose = FALSE)

# Find new biomarkers
biomarkers <- FindAllMarkers(sobj, only.pos = TRUE, verbose = FALSE)

# Save out files
## UMAP
umap <- as_tibble(
    as.data.frame(Embeddings(sobj, reduction = "umap")) %>%
        rownames_to_column("cell_id")
)
write_feather(umap, path = snakemake@output[["umap"]])

## Clusters
metadata <- as_tibble(
    sobj@meta.data %>%
        rownames_to_column("cell_id") %>%
        mutate(cluster = seurat_clusters) %>%
        select(cell_id, orig.ident, nCount_RNA, nFeature_RNA, sample, cluster)
)
write_feather(metadata, path = snakemake@output[["metadata"]])

## Biomarkers
biomarkers <- as_tibble(biomarkers) %>%
    rename(gene_symbol = gene) %>%
    left_join(fbgn2symbol) %>%
    select(FBgn, gene_symbol, cluster, everything())
write_feather(biomarkers, path = snakemake@output[["biomarkers"]])