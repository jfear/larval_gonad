library("dplyr")
library("ggplot2")
library("cowplot")
library("pheatmap")
library("RColorBrewer")
library("DESeq2")

INPUT_FILE <- "../config/sampletable.tsv"
INPUT_COUNTS <- "../../output/neox-wf/raw_counts.feather"

sample_table <- read.csv(INPUT_FILE, sep = "\t", stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(species = stringr::str_match(title, "(\\w+)_\\w+_\\w+_[\\w\\d]+")[[2]]) %>%
    mutate(tissue = stringr::str_match(title, "\\w+_(\\w+)_\\w+_[\\w\\d]+")[[2]]) %>%
    mutate(sex = stringr::str_match(title, "\\w+_\\w+_(\\w+)_[\\w\\d]+")[[2]]) %>%
    mutate(rep = stringr::str_match(title, "\\w+_\\w+_\\w+_([\\w\\d]+)")[[2]]) %>%
    select(title, species, tissue, sex, rep) %>%
    filter(species %in% c("dmel", "dpse", "dvir", "dwil", "dmoj")) %>%
    filter(tissue %in% c("WB", "GO")) %>%
    tibble::column_to_rownames("title")

cnts <- feather::read_feather(INPUT_COUNTS, columns = c("FBgn", rownames(sample_table))) %>%
    replace(is.na(.), 0) %>%
    tibble::column_to_rownames("FBgn") %>%
    filter(rowSums(.) > 1)

dds <- DESeqDataSetFromMatrix(countData = cnts, colData = sample_table, design = ~ sex + species + tissue)
dds <- estimateSizeFactors(dds)

# VST Normalization
vsd <- vst(dds, blind = FALSE)

# Sample Distances
sampleDist <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sex, vsd$species, vsd$tissue, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
)

# PCA
plotPCA(vsd, intgroup = c("sex", "species", "tissue"))

# Differential Expression
dds <- DESeq(dds)

res <- results(dds, contrast = c("species", "sex", "WB"))