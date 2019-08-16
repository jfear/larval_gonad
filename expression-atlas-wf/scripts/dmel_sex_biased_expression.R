library(DESeq2)
library(tibble)
library(dplyr)
library(ggplot2)
set.seed(42)

SAMPLE_TABLE <- snakemake@input[["sample_table"]]
COUNTS_TABLE <- snakemake@input[["counts_table"]]
GENE_ANNOTATION <- snakemake@input[["gene_annot"]]
OUTPUT_FILE <- snakemake@output[[1]]
LOG <- file(snakemake@log[[1]], open = "wt")
sink(LOG)
sink(LOG, type = "message")

OUT_DIR <- normalizePath(dirname(OUTPUT_FILE))

# Debug Settings
# SAMPLE_TABLE <- "../config/sampletable.tsv"
# COUNTS_TABLE <- "../../output/expression-atlas-wf/raw_counts.feather"
# GENE_ANNOTATION <- "../../references/gene_annotation_dmel_r6-24.feather"

# Load metadata
sampletable <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
  filter(tissue == "whole body") %>%
  filter(species == "Drosophila melanogaster") %>%
  mutate(sn = samplename) %>%
  mutate(sex = as.factor(sex)) %>%
  mutate(species = as.factor(species)) %>%
  tidyr::separate(sn, "species_abbrev", extra = "drop") %>%
  mutate(species_abbrev = as.factor(species_abbrev)) %>%
  column_to_rownames("samplename")


fbgn2symbol <- feather::read_feather(
  GENE_ANNOTATION,
  columns = c("FBgn", "gene_symbol")
)

# Load count data
raw_counts <- feather::read_feather(
  COUNTS_TABLE,
  columns = c("FBgn", row.names(sampletable))
) %>%
  filter(rowSums(select_if(., is.numeric)) != 0) %>%
  arrange(FBgn) %>%
  column_to_rownames("FBgn")


# Explority Analysis
dds <- DESeqDataSetFromMatrix(
  raw_counts[, row.names(sampletable)],
  colData = sampletable,
  design = ~sex
)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
print(sizeFactors(dds))

## Plot sample similarity
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- NULL
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

p <- pheatmap::pheatmap(sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col = colors
)

svg(file.path(OUT_DIR, "sex_similarity.svg"), width = 8, height = 8)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()


## Plot PCA
plotPCA(vsd, intgroup = c("sex", "species_abbrev")) +
  geom_text(
    aes(label = ifelse(PC2 > 0, colnames(vsd), "")),
    hjust = -0.1,
    vjust = 0.5
  ) +
  coord_fixed(ratio = 5)
ggsave(file.path(OUT_DIR, "sex_pca.svg"), width = 8, height = 8)

# Differential expression analysis
dds <- DESeq(dds, betaPrior = FALSE, fitType = "parametric")
print(resultsNames(dds))

res <- results(dds, contrast = c("sex", "male", "female"), alpha = 0.01)

res.lfc <- lfcShrink(
  dds,
  contrast = c("sex", "male", "female"),
  res = res,
  type = "normal"
)

head(res.lfc)
summary(res.lfc)

# Save results
deg_res <- as.data.frame(res.lfc) %>% rownames_to_column("FBgn")
write.table(
  deg_res,
  file = OUTPUT_FILE,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
