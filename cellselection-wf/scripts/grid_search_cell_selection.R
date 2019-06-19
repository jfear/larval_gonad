library(Seurat)
library(dplyr)
library(feather)

SAMPLE <- snakemake@wildcards[['sample']]
DATA_DIR <- dirname(snakemake@input[['mtx']])
cell_calls <- snakemake@input[['cell_calls']]
doublets <- snakemake@input[['doublets']]
OUTPUT <- snakemake@output[[1]]

# Parameters for Grid Search
check <- function(x){
    if (x == 'None'){
        return(Inf)
    }
    return(as.integer(x))
}

gene_low <- check(snakemake@wildcards[['gene_low']])
gene_high <- check(snakemake@wildcards[['gene_high']])
pct_mito <- check(snakemake@wildcards[['pct_mito']])
pct_ribo <- check(snakemake@wildcards[['pct_ribo']])

RESOLUTIONS <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2)

# Basic Seurat Clustering using 12 PCs
tenX.data <- Read10X(data.dir = DATA_DIR)
sobj <- CreateSeuratObject(
    raw.data = tenX.data,
    min.cells = 3,
    min.genes = 200,
    project = SAMPLE,
    display.progress = FALSE
)

# Get the percent of mito and rRNA genes
mito.genes <- grep(pattern = "^mt:", x = rownames(x = sobj@data), value = TRUE)
percent.mito <- Matrix::colSums(sobj@raw.data[mito.genes, ])/Matrix::colSums(sobj@raw.data)
sobj <- AddMetaData(object = sobj, metadata = percent.mito, col.name = "percent.mito")

ribo.genes <- grep(pattern = "^[^mt].*rRNA.*", x = rownames(x = sobj@data), value = TRUE)
percent.ribo <- Matrix::colSums(sobj@raw.data[ribo.genes, ])/Matrix::colSums(sobj@raw.data)
sobj <- AddMetaData(object = sobj, metadata = percent.ribo, col.name = "percent.ribo")

# Filter based on the number of genes
sobj <- FilterCells(
    object = sobj,
    subset.names = "nGene",
    low.thresholds = gene_low,
    high.thresholds = gene_high
)

# Filter based on percent mitochondria
sobj <- FilterCells(
    object = sobj,
    subset.names = "percent.mito",
    high.thresholds = pct_mito
)

# Filter based on percent rRNA
sobj <- FilterCells(
    object = sobj,
    subset.names = "percent.ribo",
    high.thresholds = pct_ribo
)

# Normalization
sobj <- NormalizeData(
    object = sobj,
    normalization.method = "LogNormalize",
    scale.factor = 1e4,
    display.progress = FALSE
)

sobj <- FindVariableGenes(
    object = sobj,
    do.plot = TRUE,
    display.progress = FALSE
)

sobj <- ScaleData(
    object = sobj,
    vars.to.regress = "nUMI",
    display.progress = FALSE
)

# Dimensionality reduction
sobj <- RunPCA(
    object = sobj, 
    pc.genes = sobj@var.genes, 
    do.print = FALSE 
)

sobj <- ProjectPCA(object = sobj, do.print = FALSE)

# Clustering 
sobj <- FindClusters(
    object = sobj, 
    reduction.type = "pca", 
    resolution = RESOLUTIONS, 
    dims.use = 1:12, 
    save.SNN = FALSE, 
    print.output = FALSE
)

# Pull out first cluster call with >= 13 clusters
metadata <- tibble::as_tibble(tibble::rownames_to_column(sobj@meta.data, 'cell_id')) %>%
    mutate(gene_low_threshold = gene_low, gene_high_threshold = gene_high, percent_mito_threshold = pct_mito, percent_ribo_threshold = pct_ribo)

for (i in 1:length(RESOLUTIONS)){
    res <- paste0('res.', RESOLUTIONS[[i]])
    lvls <- metadata[res] %>% distinct() %>% pull()
    if (length(lvls) >= 13){
        df <- metadata %>% 
            select("cell_id", "nGene", "nUMI", percent_mito = "percent.mito", percent_ribo = "percent.ribo", "gene_low_threshold", "gene_high_threshold", "percent_mito_threshold", "percent_ribo_threshold", cluster = res)
        break
    } else if (i == length(RESOLUTIONS)){
        df <- metadata %>% 
            select("cell_id", "nGene", "nUMI", percent_mito = "percent.mito", percent_ribo = "percent.ribo", "gene_low_threshold", "gene_high_threshold", "percent_mito_threshold", "percent_ribo_threshold", cluster = res)
    }
}

# Write out table with all results
write_feather(df, OUTPUT)