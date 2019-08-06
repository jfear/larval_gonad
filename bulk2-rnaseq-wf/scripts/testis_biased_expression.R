library(DESeq2)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
set.seed(42)

# Main Folders
PROJ_DIR <- "/home/fearjm/Projects/larval_gonad"
CONFIG_DIR <- file.path(PROJ_DIR, "bulk2-rnaseq-wf/config")
DATA_DIR <- file.path(PROJ_DIR, "output/bulk2-rnaseq-wf")
REF_DIR <- file.path(PROJ_DIR, "references")

# Main Files
SAMPLE_TABLE <- file.path(CONFIG_DIR, "sampletable.tsv")
COUNTS_TABLE <- file.path(DATA_DIR, "rnaseq_aggregation/gene_level_counts.tsv")
INTERGENIC_COUNTS_TABLE <- file.path(DATA_DIR, "rnaseq_aggregation/intergenic_counts.tsv")
GENE_ANNOTATION <- file.path(REF_DIR, "gene_annotation_dmel_r6-24.feather")
ERCC_ANNOTATION <- file.path(REF_DIR, "ercc_annotation.tsv")


# Load metadata
sampletable <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
    tibble::column_to_rownames("samplename")

fbgn2symbol <- feather::read_feather(
    GENE_ANNOTATION,
    columns = c("FBgn", "gene_symbol")
)

ercc_conc <- read.csv(ERCC_ANNOTATION, stringsAsFactors = FALSE, sep = "\t") %>%
    mutate(subpool = Subpool_in_pool_78, Aconc = X78A_nmol_per_ul, Bconc = X78B_nmol_per_ul) %>%
    select(ercc_id, subpool, Aconc, Bconc)

# Load raw count data
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

combined_counts <- rbind(ercc_counts, raw_counts, intergenic_counts) %>% tibble::column_to_rownames("Geneid")

# Normalize combined counts
## We use all counts (ercc, gene, intergenic) to estimate size factors and normalize
dds <- DESeqDataSetFromMatrix(combined_counts[, row.names(sampletable)], sampletable, design = ~group)
dds <- estimateSizeFactors(dds)
print(sizeFactors(dds))

norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame()

# Calculate Intergenic Cutoff
## The low expression threshold is defined in:
## > Lin, Yanzhu, Kseniya Golovnina, Zhen-Xia Chen, Hang Noh Lee, Yazmin L.
## > Serrano Negron, Hina Sultana, Brian Oliver, and Susan T. Harbison. 2016.
## > “Comparison of Normalization and Differential Expression Analyses Using
## > RNA-Seq Data from 726 Individual Drosophila Melanogaster.” BMC Genomics 17
## > (1): 28.
##
## We pull out intergenic reads and find the log value at the 95th percentile. I
## am not sure why they do all of this log conversion. If you do the same thing 
## without the log you get the about same answer.
LET <- norm_counts %>%
    filter(stringr::str_detect(row.names(.), "^inter")) %>%
    purrr::map(function(x) quantile(log2(x + 1), 0.95))

LET <- quantile(as.matrix(log2(norm_counts[grep("inter", row.names(norm_counts)),] + 1)), 0.95)

## We then get convert back out of log2 space.
oLET <- lapply(LET, function(x) 2^x - 1)
print(as.data.frame(oLET))
oLET <- 2^LET - 1

## Now we have a low expression threshold `oLET`. We will set any values less
## than this to 0.

# Set all counts below oLET to 0
zero_data <- function(df_norm, df_raw, oLET, pattern){
    # Create as mask for values <= oLET
    norm <- df_norm %>%
        tibble::rownames_to_column("Geneid") %>%
        filter(stringr::str_detect(Geneid, pattern)) %>%
        arrange(Geneid) %>%
        tibble::column_to_rownames("Geneid")

    mask <- norm <= oLET

    # Take raw counts table and set values <= oLET to 0
    counts <- df_raw %>% 
        arrange(Geneid) %>% 
        tibble::column_to_rownames("Geneid")

    counts[mask] = 0

    # Drop rows that are all 0
    return(counts[rowSums(counts) != 0, ])
}

fbgn_zcounts <- zero_data(norm_counts, raw_counts, oLET, "^FBgn")
ercc_zcounts <- zero_data(norm_counts, ercc_counts, oLET, "^ERCC")

# ERCC Analysis
ercc_zcounts_melted <- ercc_zcounts %>%
    tibble::rownames_to_column("ercc_id") %>%
    tidyr::gather("samplename", "count", -ercc_id)

sample2pool <- sampletable %>% 
    tibble::rownames_to_column("samplename") %>% 
    mutate(pool = ercc) %>% 
    select(samplename, pool)

ercc_stacked <- ercc_zcounts_melted %>%
    left_join(sample2pool) %>%
    group_by(ercc_id, pool) %>%
    summarize(total_count = sum(count)) %>%
    tidyr::spread(pool, total_count) %>%
    filter(A > 0 & B > 0)

# Plot counts vs conc
conc_data <- ercc_stacked %>% 
    left_join(ercc_conc)

## Pool A
ggplot(conc_data, aes(x=log10(Aconc), y=log10(A))) + 
    geom_point() + 
    geom_smooth(se = FALSE, method="lm", na.rm = TRUE, color="black") +
    geom_text(aes(label=ercc_id), hjust=0, vjust=0) + 
    xlab("log10 Concentration") +
    ylab("log10 Read Counts") + 
    ggtitle("Pool A")

## Pool B
ggplot(conc_data, aes(x=log10(Bconc), y=log10(B))) + 
    geom_point() + 
    geom_smooth(se = FALSE, method="lm", na.rm = TRUE, color="black") +
    geom_text(aes(label=ercc_id), hjust=0, vjust=0) + 
    xlab("log10 Concentration") +
    ylab("log10 Read Counts") + 
    ggtitle("Pool B")

# Plot ratios of pool A / B
ercc_ratio <- ercc_stacked %>%
    mutate(
        baseMean = mean(c(A, B)),
        logFC = log2(A / B)
    ) %>%
    select(ercc_id, baseMean, logFC) %>%
    left_join(ercc_conc)

ggplot(ercc_ratio, aes(x=log2(baseMean), y=logFC, color=subpool)) + 
    geom_point() + 
    ylim(-1, 1) + 
    xlim(4, 18) + 
    geom_smooth(se = FALSE, method="loess", na.rm = TRUE) +
    geom_text(aes(label=ercc_id), hjust=0, vjust=0)


# Differential expression analysis
dds <- DESeqDataSetFromMatrix(fbgn_zcounts[, row.names(sampletable)], colData=sampletable, design= ~group)
dds <- estimateSizeFactors(dds)
print(sizeFactors(dds))

vsd <- vst(dds, blind=FALSE)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap::pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

plotPCA(vsd, intgroup = "group") + geom_text(aes(label=colnames(vsd)), hjust=-0.1, vjust=0.5)

dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds, contrast=c("group", "TCP", "OCP"), alpha=0.01)
res.lfc <- lfcShrink(dds, contrast=c("group", "TCP", "OCP"), res=res)
summary(res)
summary(res.lfc)