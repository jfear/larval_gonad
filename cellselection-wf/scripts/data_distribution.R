library(Seurat)
library(dplyr)
library(ggplot2)
library(feather)
source("scripts/common.R")

SAMPLE <- snakemake@wildcards[['sample']]
DATA_DIR <- dirname(snakemake@input[['mtx']])
CELL_CALLS <- snakemake@input[['cell_calls']]
DOUBLETS <- snakemake@input[['doublets']]
GENE_ANNOT <- snakemake@input[['gene_annotation']]
VLN <- snakemake@output[['vln']]
SCATTER <- snakemake@output[['scatter']]
REP <- gsub('testis', 'rep', SAMPLE)

# Get list of good cell IDs
cells <- read_feather(CELL_CALLS, columns=c('cell_id', 'is_cell')) %>% 
    filter(.$is_cell) %>%
    pull(cell_id)

dbls <- read.table(DOUBLETS, stringsAsFactors = FALSE, header = FALSE, col.names = c('cell_id')) %>% 
    pull(cell_id)

good_cells <- cells[!cells %in% dbls]

# Gene annotation
fbgn2symbol <- read_feather(GENE_ANNOT, columns=c('FBgn', 'gene_symbol'))

# Read and filter good cells from 10X
tenX.data <- load_and_filter_10x(DATA_DIR, good_cells, fbgn2symbol, REP)

# Build basic seurat object
sobj <- CreateSeuratObject(
    raw.data = tenX.data,
    min.cells = 3,
    project = SAMPLE,
    display.progress = FALSE
)

# Add on mitochondrial and rRNA annotations for plotting
mito.genes <- grep(pattern = "^mt:", x = rownames(x = sobj@data), value = TRUE)
percent.mito <- Matrix::colSums(sobj@raw.data[mito.genes, ])/Matrix::colSums(sobj@raw.data)
sobj <- AddMetaData(object = sobj, metadata = percent.mito, col.name = "percent.mito")

ribo.genes <- grep(pattern = "^[^mt].*rRNA.*", x = rownames(x = sobj@data), value = TRUE)
percent.ribo <- Matrix::colSums(sobj@raw.data[ribo.genes, ])/Matrix::colSums(sobj@raw.data)
sobj <- AddMetaData(object = sobj, metadata = percent.ribo, col.name = "percent.ribo")

# Distribution plot
plt <- VlnPlot(
    object = sobj,
    features.plot = c("nGene", "nUMI", "percent.mito", "percent.ribo"),
    nCol = 4,
    point.size.use = .05
)
ggsave(VLN, plot = plt, width = 8, height = 6)

# Scatter plot
.cor <- round(cor(sobj@meta.data$nGene, sobj@meta.data$nUMI), 2)
p1 <- ggplot(sobj@meta.data, aes(x=nUMI, y=nGene)) + geom_point() + scale_x_log10() + scale_y_log10() + xlab("UMI") + ylab("Genes") + ggtitle(paste0(SAMPLE, ': (', .cor, ')'))

.cor <- round(cor(sobj@meta.data$percent.mito, sobj@meta.data$nUMI), 2)
p2 <- ggplot(sobj@meta.data, aes(x=nUMI, y=percent.mito)) + geom_point() + scale_x_log10() + xlab("UMI") + ylab("Percent Mitochondrial") + ggtitle(paste0(SAMPLE, ': (', .cor, ')'))

.cor <- round(cor(sobj@meta.data$percent.ribo, sobj@meta.data$nUMI), 2)
p3 <- ggplot(sobj@meta.data, aes(x=nUMI, y=percent.ribo)) + geom_point() + scale_x_log10() + xlab("UMI") + ylab("Percent Ribosomal RNA") + ggtitle(paste0(SAMPLE, ': (', .cor, ')'))
plt <- gridExtra::grid.arrange(p1, p2, p3)
ggsave(SCATTER, plot = plt, width = 8, height = 12)
