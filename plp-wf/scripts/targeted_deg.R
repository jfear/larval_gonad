library(Seurat)
library(tibble)
library(dplyr)
library(feather)
set.seed(42)

# Load Gene Metadata
fbgn2symbol <- read_feather(
    "../../references/gene_annotation_dmel_r6-26.feather",
    columns = c("FBgn", "gene_symbol")
)

# Load Combined Data
load("../../output/seurat3-cluster-wf/combined_n3.Robj")

# Run PCA
pca <- combined@reductions$pca@cell.embeddings[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("cell_id")
write.table(pca, "../../output/plp-wf/pca.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Trial List
trial_fbgns <- read.table(
    "../../data/external/galletta/trial_list.txt",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
) %>% pull(FBGN)

trial_symbols <- fbgn2symbol %>%
    filter(FBgn %in% trial_fbgns) %>%
    pull(gene_symbol)

# Pull out trial genes
sobj <- subset(x = combined, features = trial_symbols)

# G vs All
deg <- FindMarkers(sobj, "6", c("4", "0", "2"), logfc.threshold = 0, min.pct = 0.01)
deg_clean <- as_tibble(deg %>% rownames_to_column("gene_symbol")) %>%
    left_join(fbgn2symbol) %>%
    select(FBgn, gene_symbol, p_val, p_val_adj, everything()) %>%
    arrange(desc(avg_logFC))

write.table(deg_clean, file = "../../output/plp-wf/trial_list_deg_gonia_vs_spermatocytes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# G vs EPS
deg <- FindMarkers(sobj, "6", "4", logfc.threshold = 0, min.pct = 0.01)
deg_clean <- as_tibble(deg %>% rownames_to_column("gene_symbol")) %>%
    left_join(fbgn2symbol) %>%
    select(FBgn, gene_symbol, p_val, p_val_adj, everything()) %>%
    arrange(desc(avg_logFC))

write.table(deg_clean, file = "../../output/plp-wf/trial_list_deg_gonia_vs_eps.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# G vs MPS and LPS
deg <- FindMarkers(sobj, "6", c("0", "2"), logfc.threshold = 0, min.pct = 0.01)
deg_clean <- as_tibble(deg %>% rownames_to_column("gene_symbol")) %>%
    left_join(fbgn2symbol) %>%
    select(FBgn, gene_symbol, p_val, p_val_adj, everything()) %>%
    arrange(desc(avg_logFC))

write.table(deg_clean, file = "../../output/plp-wf/trial_list_deg_gonia_vs_mid_and_late.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# G vs MPS
deg <- FindMarkers(sobj, "6", "0", logfc.threshold = 0, min.pct = 0.01)
deg_clean <- as_tibble(deg %>% rownames_to_column("gene_symbol")) %>%
    left_join(fbgn2symbol) %>%
    select(FBgn, gene_symbol, p_val, p_val_adj, everything()) %>%
    arrange(desc(avg_logFC))

write.table(deg_clean, file = "../../output/plp-wf/trial_list_deg_gonia_vs_mid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# G vs LPS
deg <- FindMarkers(sobj, "6", "2", logfc.threshold = 0, min.pct = 0.01)
deg_clean <- as_tibble(deg %>% rownames_to_column("gene_symbol")) %>%
    left_join(fbgn2symbol) %>%
    select(FBgn, gene_symbol, p_val, p_val_adj, everything()) %>%
    arrange(desc(avg_logFC))

write.table(deg_clean, file = "../../output/plp-wf/trial_list_deg_gonia_vs_late.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
