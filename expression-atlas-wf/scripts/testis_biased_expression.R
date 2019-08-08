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
# SAMPLE_TABLE <- "../config/sampletable.tsv"
# COUNTS_TABLE <- "../../output/expression-atlas-wf/raw_counts.feather"
# GENE_ANNOTATION <- "../../references/gene_annotation_dmel_r6-24.feather"

# Load metadata
sampletable <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
  filter(tissue == "gonad") %>%
  filter(species == "Drosophila melanogaster") %>%
  mutate(sn = samplename) %>%
  tidyr::separate(sn, "species_abbrev", extra = "drop") %>%
  mutate(sex = as.factor(sex)) %>%
  mutate(species_abbrev = as.factor(species_abbrev)) %>%
  column_to_rownames("samplename")


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
dds <- DESeqDataSetFromMatrix(raw_counts[, row.names(sampletable)], colData = sampletable, design = ~ sex*species_abbrev)
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
dds <- DESeq(dds, betaPrior = FALSE, fitType = "parametric")

res_sex <- results(dds, contrast = c("sex", "male", "female"), alpha = 0.01)
head(res_sex)
summary(res_sex)

res_species <- results(dds, contrast = c("species_abbrev", "w1118", "orgR"), alpha = 0.01)
head(res_species)
summary(res_species)

res_interact <- results(dds, name = "sexmale.species_abbrevw1118", alpha = 0.01)
head(res_interact)
summary(res_interact)

res.lfc <- lfcShrink(dds, contrast = c("sex", "male", "female"), res = res, type = "normal")
head(res.lfc)
summary(res.lfc)

# Save results
deg_res <- as.data.frame(res.lfc) %>% rownames_to_column("FBgn")
write.table(deg_res, file = OUTPUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)