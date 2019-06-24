
library(Seurat)
library(dplyr)
library(feather)
set.seed(42)

INPUT_ROBJ <- gsub(".html", ".Robj", snakemake@input[["html"]])
GENE_ANNOTATION <- snakemake@input[["gene_annotation"]]

TSV <- snakemake@output[['tsv']]
FEATHER <- snakemake@output[['feather']]

# Load the data
load(INPUT_ROBJ)

fbgn2symbol <- read_feather(
    GENE_ANNOTATION,
    columns = c("FBgn", "gene_symbol")
)

# Seurat suggests to use the RNA assays for DEG
DefaultAssay(combined) <- "RNA"

# Call biomarkers
biomarkers <- FindAllMarkers(combined, only.pos = TRUE, verbose = FALSE)

# Merge on FBgn
biomarkers <- tibble::as_tibble(biomarkers) %>%
    rename(gene_symbol = gene) %>%
    left_join(fbgn2symbol) %>%
    select(FBgn, gene_symbol, cluster, everything())

# Save outputs
write.table(
    biomarkers,
    file = TSV,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

write_feather(
    biomarkers,
    path = FEATHER
)