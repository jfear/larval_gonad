library(Seurat)
library(dplyr)
library(ggplot)

SAMPLE <- snakemake@wildcards[['sample']]
DATA_DIR <- dirname(snakemake@input[['mtx']])
cell_calls <- snakemake@input[['cell_calls']]
doublets <- snakemake@input[['doublets']]
VLN <- snakemake@output[['vln']]
SCATTER <- snakemake@output[['scatter']]

# Basic Seurat Clustering using 12 PCs
tenX.data <- Read10X(data.dir = DATA_DIR)
sobj <- CreateSeuratObject(
    raw.data = tenX.data,
    min.cells = 3,
    min.genes = 200,
    project = SAMPLE,
    display.progress = FALSE
)
