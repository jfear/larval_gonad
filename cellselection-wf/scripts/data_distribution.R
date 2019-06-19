library(Seurat)
library(dplyr)
library(ggplot2)
library(feather)

# SAMPLE <- snakemake@wildcards[['sample']]
# DATA_DIR <- dirname(snakemake@input[['mtx']])
# cell_calls <- snakemake@input[['cell_calls']]
# doublets <- snakemake@input[['doublets']]
# VLN <- snakemake@output[['vln']]
# SCATTER <- snakemake@output[['scatter']]

SAMPLE <- 'testis1'
DATA_DIR <- '../../output/cellranger-wf/testis1/outs/raw_gene_bc_matrices/dmelr6-24'
cell_calls <- '../../output/cellselection-wf/dropletutils/testis1_cell_calls.feather'
doublets <- '../../output/cellselection-wf/testis1_scrublet_dublets.txt'
GENE_ANNOT <- '../../references/gene_annotation_dmel_r6-24.feather'
VLN <- '../../output/cellselection-wf/testis1_distribtuion_vln.png'
SCATTER <- '../../output/cellselection-wf/testis1_distribtuion_scatter.png'

REP <- gsub('testis', 'rep', SAMPLE)

# Get list of good cell IDs
cells <- read_feather(cell_calls, columns=c('cell_id', 'is_cell')) %>% 
    filter(.$is_cell) %>%
    pull(cell_id)

dbls <- read.table(doublets, stringsAsFactors = FALSE, header = FALSE, col.names = c('cell_id')) %>% 
    pull(cell_id)

good_cells <- cells[!cells %in% dbls]

# Load 10x data
tenX.data <- Read10X(data.dir = DATA_DIR)

# Keep good cells
colnames(tenX.data) <- as.vector(sapply(colnames(tenX.data), function(x) paste(REP, x, sep='_')))
tenX.data <- tenX.data[, good_cells]

# Rename to gene symbol
fbgn2symbol <- read_feather(GENE_ANNOT, columns=c('FBgn', 'gene_symbol'))
row.names(tenX.data) <- left_join(tibble(FBgn = row.names(tenX.data)), fbgn2symbol, by='FBgn') %>% pull(gene_symbol)

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
