library(Seurat)
library(dplyr)
library(feather)
library(tibble)
set.seed(42)

HTML <- snakemake@input[["html"]]
if (grepl("combined", HTML)) {
    COMBINED <- TRUE
    INPUT_ROBJ <- gsub(".html", ".Robj", HTML)
} else {
    COMBINED <- FALSE
    INPUT_ROBJ <- gsub("_individual.html", ".Robj", HTML)
}

GENE_ANNOTATION <- snakemake@input[["gene_annotation"]]

NORM <- snakemake@output[["norm"]]
METADATA <- snakemake@output[["metadata"]]
PCA_EMBED <- snakemake@output[["pca_embed"]]
PCA_LOAD <- snakemake@output[["pca_load"]]
UMAP <- snakemake@output[["umap"]]

RESOLUTION <- snakemake@params[[1]]

# Load the data
load(INPUT_ROBJ)
if (COMBINED) {
    sobj <- combined
    DefaultAssay(sobj) <- "RNA"
}

fbgn2symbol <- read_feather(
    GENE_ANNOTATION,
    columns = c("FBgn", "gene_symbol")
)

# Normalized data
df <- as_tibble(
    as.data.frame(GetAssayData(sobj)) %>%
        rownames_to_column("gene_symbol")
) %>%
    left_join(fbgn2symbol, by = "gene_symbol") %>%
    select(FBgn, gene_symbol, everything())

write_feather(df, path = NORM)

# Meatadata
df <- as_tibble(
    sobj@meta.data %>%
        rownames_to_column("cell_id")
) %>%
    select(cell_id, nCount_RNA, nFeature_RNA, starts_with("int"))

df$cluster <- df[, RESOLUTION] %>% pull(.)

write_feather(df, path = METADATA)

# PCA Cell Embedings
df <- as_tibble(
    as.data.frame(Embeddings(sobj, reduction = "pca")) %>%
        rownames_to_column("cell_id")
)

write_feather(df, path = PCA_EMBED)

# PCA Gene Loadings
df <- as_tibble(
    as.data.frame(Loadings(sobj, reduction = "pca")) %>%
        rownames_to_column("gene_symbol") %>%
        left_join(fbgn2symbol, by = "gene_symbol") %>%
        select(FBgn, gene_symbol, everything())
)

write_feather(df, path = PCA_LOAD)

# UMAP
df <- as_tibble(
    as.data.frame(Embeddings(sobj, reduction = "umap")) %>%
        rownames_to_column("cell_id")
)
write_feather(df, path = UMAP)