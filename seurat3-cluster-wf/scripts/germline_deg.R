library(Seurat)
library(tibble)
library(dplyr)
library(feather)
set.seed(42)

# Load Gene Metadata
fbgn2symbol <- feather::read_feather(
    snakemake@input[["gene_annotation"]],
    columns = c("FBgn", "gene_symbol")
)

# Load Combined Data
robj <- gsub(".html", ".Robj", snakemake@input[["html"]])
load(robj)

# Helper function
find_markers <- function(obj, ident.1, ident.2, alpha = 1) {
    deg <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2)

    deg_clean <- as_tibble(deg %>% tibble::rownames_to_column("gene_symbol")) %>%
        left_join(fbgn2symbol) %>%
        select(FBgn, gene_symbol, everything()) %>%
        filter(p_val_adj <= alpha) %>%
        arrange(desc(avg_logFC))

    return(deg_clean)
}

# GvPS
df <- find_markers(combined, "6", c("0", "2", "4"))
write_feather(df, path = snakemake@output[["gvps"]])

# GvEPS
df <- find_markers(combined, "6", "4")
write_feather(df, path = snakemake@output[["gveps"]])

# GvMLPS
df <- find_markers(combined, "6", c("0", "2"))
write_feather(df, path = snakemake@output[["gvmlps"]])

# GvMPS
df <- find_markers(combined, "6", "0")
write_feather(df, path = snakemake@output[["gvmps"]])

# GvLPS
df <- find_markers(combined, "6", "2")
write_feather(df, path = snakemake@output[["gvlps"]])

# EPSvMLPS
df <- find_markers(combined, "4", c("0", "2"))
write_feather(df, path = snakemake@output[["epsvmlps"]])

# EPSvMPS
df <- find_markers(combined, "4", "0")
write_feather(df, path = snakemake@output[["epsvmps"]])

# EPSvLPS
df <- find_markers(combined, "4", "2")
write_feather(df, path = snakemake@output[["epsvlps"]])

# MPSvLPS
df <- find_markers(combined, "0", "2")
write_feather(df, path = snakemake@output[["mpsvlps"]])
