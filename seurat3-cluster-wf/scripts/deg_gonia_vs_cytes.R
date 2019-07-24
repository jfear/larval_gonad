library(Seurat)
library(dplyr)
set.seed(42)

ROBJ <- gsub(".html", ".Robj", snakemake@input[["html"]])
GENE_METADATA <- snakemake@input[["gene_annotation"]]

GVC <- list(feather = snakemake@output[["gvcf"]], tsv = snakemake@output[["gvct"]])
GVE <- list(feather = snakemake@output[["gvef"]], tsv = snakemake@output[["gvet"]])
EVP <- list(feather = snakemake@output[["evpf"]], tsv = snakemake@output[["evpt"]])

# Debug Settings
# ROBJ <- "../../output/seurat3-cluster-wf/combined_n3.Robj"
# GENE_METADATA <- "../../references/gene_annotation_dmel_r6-24.feather"
# GVC <- list(feather = "/tmp/test.feather", tsv = "/tmp/test.tsv")

# Load Gene Metadata
fbgn2symbol <- feather::read_feather(GENE_METADATA, columns = c("FBgn", "gene_symbol"))

# Load Combined Data
load(ROBJ)

# Helper function
find_markers <- function(obj, ident.1, ident.2, alpha = 1) {
    deg <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2)

    deg_clean <- tibble::as_tibble(deg %>% tibble::rownames_to_column("gene_symbol")) %>%
        left_join(fbgn2symbol) %>%
        select(FBgn, gene_symbol, everything()) %>%
        filter(p_val_adj <= alpha) %>%
        arrange(desc(avg_logFC))

    return(deg_clean)
}

save_deg <- function(df, output_names) {
    write.table(
        df,
        file = output_names[["tsv"]],
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
    )

    feather::write_feather(
        df,
        path = output_names[["feather"]]
    )
}

gonia_vs_cytes <- find_markers(combined, "9", c("0", "2", "5", "6"))
save_deg(gonia_vs_cytes, GVC)

gonia_vs_eps <- find_markers(combined, "9", "5")
save_deg(gonia_vs_eps, GVE)

eps_vs_ps <- find_markers(combined, "5", c("0", "2", "6"))
save_deg(eps_vs_ps, EVP)
