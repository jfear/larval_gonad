library(dplyr)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")
library(DESeq2)


# Import raw counts
raw_counts <- feather::read_feather('../../output/neox-wf/raw_counts.feather') %>% 
    tibble::column_to_rownames('FBgn') %>% 
    replace(., is.na(.), 0)

# Import sample table and extract useful design info
sample_table <- read.csv('../config/sampletable.tsv', sep='\t', stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(species = strsplit(title, split = "_")[[1]][[1]]) %>%
    mutate(tissue = strsplit(title, split = "_")[[1]][[2]]) %>%
    mutate(sex = strsplit(title, split = "_")[[1]][[3]]) %>%
    mutate(rep = strsplit(title, split = "_")[[1]][[4]]) %>% 
    tibble::column_to_rownames('title')

# Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = sample_table, design = ~ species + tissue + sex)

# Remove rows of 0
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# Different Normalizations 
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

# Plot normalizations
dds <- estimateSizeFactors(dds)

df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
        mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)  




