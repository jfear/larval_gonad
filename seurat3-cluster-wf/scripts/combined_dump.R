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
    select(cell_id, nCount_RNA, nFeature_RNA, RESOLUTION) %>%
    rename(cluster = RESOLUTION)

write_feather(df, path = METADATA)