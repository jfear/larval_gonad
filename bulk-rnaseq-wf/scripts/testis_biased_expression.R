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
# CONFIG_DIR <- file.path(PROJ_DIR, "bulk-rnaseq-wf/config")
# DATA_DIR <- file.path(PROJ_DIR, "output/bulk-rnaseq-wf")
# REF_DIR <- file.path(PROJ_DIR, "references")
# SAMPLE_TABLE <- file.path(CONFIG_DIR, "sampletable.tsv")
# COUNTS_TABLE <- file.path(DATA_DIR, "rnaseq_aggregation/gene_level_counts.tsv")
# INTERGENIC_COUNTS_TABLE <- file.path(DATA_DIR, "rnaseq_aggregation/intergenic_counts.tsv")
# GENE_ANNOTATION <- file.path(REF_DIR, "gene_annotation_dmel_r6-26.feather")
# ERCC_ANNOTATION <- file.path(REF_DIR, "ercc_annotation.tsv")

# Load metadata
sampletable <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE, sep = "\t", comment.char = "#") %>%
  column_to_rownames("samplename")
sampletable$group <- as.factor(sampletable$group)

fbgn2symbol <- feather::read_feather(
  GENE_ANNOTATION,
  columns = c("FBgn", "gene_symbol")
)

ercc_conc <- read.csv(ERCC_ANNOTATION, stringsAsFactors = FALSE, sep = "\t") %>%
  mutate(subpool = Subpool_in_pool_78, Aconc = X78A_nmol_per_ul, Bconc = X78B_nmol_per_ul) %>%
  select(ercc_id, subpool, Aconc, Bconc)

# Load count data
ercc_counts <- read.csv(COUNTS_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
  filter(stringr::str_detect(Geneid, "^ERCC")) %>%
  filter(rowSums(select_if(., is.numeric)) != 0) %>%
  arrange(Geneid)

raw_counts <- read.csv(COUNTS_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
  filter(stringr::str_detect(Geneid, "^FBgn")) %>%
  filter(rowSums(select_if(., is.numeric)) != 0) %>%
  arrange(Geneid)

intergenic_counts <- read.csv(INTERGENIC_COUNTS_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
  filter(rowSums(select_if(., is.numeric)) != 0) %>%
  arrange(Geneid)

# Normalize combined counts (ercc, gene, intergenic)
combined_counts <- rbind(ercc_counts, raw_counts, intergenic_counts) %>% column_to_rownames("Geneid")

dds <- DESeqDataSetFromMatrix(combined_counts[, row.names(sampletable)], sampletable, design = ~group)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame()

print(sizeFactors(dds))

# Calculate a single 95th quantile intergenic expression cutoff according to:
# > Lin, Yanzhu, Kseniya Golovnina, Zhen-Xia Chen, Hang Noh Lee, Yazmin L.
# > Serrano Negron, Hina Sultana, Brian Oliver, and Susan T. Harbison. 2016.
# > “Comparison of Normalization and Differential Expression Analyses Using
# > RNA-Seq Data from 726 Individual Drosophila Melanogaster.” BMC Genomics 17
# > (1): 28.

LET <- quantile(as.matrix(log2(norm_counts[grep("inter", row.names(norm_counts)), ] + 1)), 0.95)

# Convert from log2 space back to normalized count space.
oLET <- 2^LET - 1
print(paste("LET: ", LET))
print(paste("oLET: ", oLET))

# Calculate a 95th quantile intergenic expression cutoff for each sample.
# NOTE: I use these values as our cutoff
LET <- norm_counts %>%
  filter(stringr::str_detect(row.names(.), "^inter")) %>%
  purrr::map(function(x) quantile(log2(x + 1), 0.95))

oLET <- lapply(LET, function(x) 2^x - 1)

let_thresh <- rbind(as.data.frame(LET), as.data.frame(oLET))
row.names(let_thresh) <- c("LET", "oLET")
print(let_thresh)

# Set reads to zero that fall below oLET.
# NOTE: oLET is for normalized reads and I want to use raw reads. I need to
# find where to zero on the normalized data and then set these to zero in the
# raw data.

## ERCCS
mask <- norm_counts[grep("ERCC", row.names(norm_counts)), ] <= oLET
print("ERCCs Below Low Expression Threshold:")
print(table(mask))

df <- ercc_counts %>% column_to_rownames("Geneid")
df.sorted <- df[match(row.names(mask), row.names(df)), match(colnames(mask), colnames(df))]
df.sorted[mask] <- 0
ercc_zcounts <- df.sorted[rowSums(df.sorted) > 0, ]  # remove rows of all 0's

## FBGNS
mask <- norm_counts[grep("FBgn", row.names(norm_counts)), ] <= oLET
print("FBgns Below Low Expression Threshold:")
print(table(mask))

df <- raw_counts %>% column_to_rownames("Geneid")
df.sorted <- df[match(row.names(mask), row.names(df)), match(colnames(mask), colnames(df))]
df.sorted[mask] <- 0
fbgn_zcounts <- df.sorted[rowSums(df.sorted) > 0, ]  # remove rows of all 0's

# ERCC Analysis
ercc_zcounts_melted <- ercc_zcounts %>%
  rownames_to_column("ercc_id") %>%
  tidyr::gather("samplename", "count", -ercc_id)

sample2pool <- sampletable %>%
  rownames_to_column("samplename") %>%
  mutate(pool = ercc) %>%
  select(samplename, pool)

ercc_stacked <- ercc_zcounts_melted %>%
  left_join(sample2pool, by = "samplename") %>%
  group_by(ercc_id, pool) %>%
  summarize(total_count = sum(count)) %>%
  tidyr::spread(pool, total_count) %>%
  filter(A > 0 & B > 0)

# Plot counts vs conc
conc_data <- ercc_stacked %>%
  left_join(ercc_conc, by = "ercc_id")

## Pool A
p1 <- ggplot(conc_data, aes(x = log10(Aconc), y = log10(A))) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", na.rm = TRUE, color = "black") +
  geom_text(aes(label = ercc_id), hjust = 0, vjust = 0) +
  xlab("log10 Concentration") +
  ylab("log10 Read Counts") +
  xlim(-10, 0) +
  ylim(.5, 5) +
  ggtitle("Pool A")

## Pool B
p2 <- ggplot(conc_data, aes(x = log10(Bconc), y = log10(B))) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", na.rm = TRUE, color = "black") +
  geom_text(aes(label = ercc_id), hjust = 0, vjust = 0) +
  xlab("log10 Concentration") +
  ylab("log10 Read Counts") +
  xlim(-10, 0) +
  ylim(.5, 5) +
  ggtitle("Pool B")

svg(file.path(OUT_DIR, "ercc_counts_vs_conc.svg"), width=8, height=8 * .8)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

# Plot ratios of pool A / B
ercc_ratio <- ercc_stacked %>%
  mutate(baseMean = mean(c(A, B)), logFC = log2(A / B)) %>%
  select(ercc_id, baseMean, logFC) %>%
  left_join(ercc_conc, by = "ercc_id")

ggplot(ercc_ratio, aes(x = log2(baseMean), y = logFC, color = subpool)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "loess", na.rm = TRUE) +
  geom_text(aes(label = ercc_id), hjust = 0, vjust = 0)

ggsave(file.path(OUT_DIR, "ercc_pool_A_vs_B.svg"), width=8, height=8 * .8)

# Explority Analysis
dds <- DESeqDataSetFromMatrix(fbgn_zcounts[, row.names(sampletable)], colData = sampletable, design = ~group)
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
svg(file.path(OUT_DIR, "bulk_sample_similarity.svg"), width=8, height=8)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()


## Plot PCA
plotPCA(vsd, intgroup = "group") +
  geom_text(aes(label = ifelse(PC2 > 0, colnames(vsd), "")), hjust = -0.1, vjust = 0.5) +
  coord_fixed(ratio = 5)
ggsave(file.path(OUT_DIR, "bulk_pca.svg"), width=8, height=8)

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
