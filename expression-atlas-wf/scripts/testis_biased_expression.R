library(DESeq2)
library(tibble)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
set.seed(42)

SAMPLE_TABLE <- snakemake@input[["sample_table"]]
COUNTS_TABLE <- snakemake@input[["counts_table"]]
INTERGENIC_COUNTS_TABLE <- snakemake@input[["intergenic_counts_table"]]
GENE_ANNOTATION <- snakemake@input[["gene_annot"]]
ERCC_ANNOTATION <- snakemake@input[["ercc_annot"]]
OUTPUT_FILE <- snakemake@output[[1]]
LOG <- file(snakemake@log[[1]], open="wt")
sink(LOG)
sink(LOG, type="message")

OUT_DIR <- normalizePath(dirname(OUTPUT_FILE))

# Debug Settings
# PROJ_DIR <- "/home/fearjm/Projects/larval_gonad"
# CONFIG_DIR <- file.path(PROJ_DIR, "expression-atlas-wf/config")
# DATA_DIR <- file.path(PROJ_DIR, "output/expression-atlas-wf")
# REF_DIR <- file.path(PROJ_DIR, "references")
# SAMPLE_TABLE <- file.path(CONFIG_DIR, "sampletable.tsv")
# COUNTS_TABLE <- file.path(DATA_DIR, "rnaseq_aggregation/gene_level_counts.tsv")
# GENE_ANNOTATION <- file.path(REF_DIR, "gene_annotation_dmel_r6-24.feather")

# Load metadata
sampletable <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
  filter(tissue = "gonad") %>%
  column_to_rownames("samplename")
sampletable$sex <- as.factor(sampletable$sex)

fbgn2symbol <- feather::read_feather(
  GENE_ANNOTATION,
  columns = c("FBgn", "gene_symbol")
)

# Load count data
raw_counts <- feather::read_feather(COUNTS_TABLE, columns = c("FBgn", row.names(sampletable))) %>%
  filter(rowSums(select_if(., is.numeric)) != 0) %>%
  arrange(FBgn) %>%
  column_to_rownames("FBgn")


# Explority Analysis
dds <- DESeqDataSetFromMatrix(raw_counts[, row.names(sampletable)], colData = sampletable, design = ~sex)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
print(sizeFactors(dds))

## Plot sample similarity
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

p <- pheatmap::pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)
svg(file.path(OUT_DIR, "gonad_similarity.svg"), width=8, height=8)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()


## Plot PCA
plotPCA(vsd, intgroup = "group") +
  geom_text(aes(label = ifelse(PC2 > 0, colnames(vsd), "")), hjust = -0.1, vjust = 0.5) +
  coord_fixed(ratio = 5)
ggsave(file.path(OUT_DIR, "gonad_pca.svg"), width=8, height=8)

# Differential expression analysis
dds <- DESeq(dds, betaPrior = FALSE, fitType = "local")
res <- results(dds, contrast = c("group", "TCP", "OCP"), alpha = 0.01)
head(res)
summary(res)

res.lfc <- lfcShrink(dds, contrast = c("group", "TCP", "OCP"), res = res, type = "normal")
head(res.lfc)
summary(res.lfc)

# Save results
deg_res <- as.data.frame(res.lfc) %>% rownames_to_column("FBgn")
write.table(deg_res, file = OUTPUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)