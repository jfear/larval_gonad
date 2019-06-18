
library(Seurat)
library(dplyr)
library(yaml)


SAMPLE <- snakemake@wildcards[['sample']]
DATA_DIR <- file.path('../', dirname(snakemake@input[["data"]]))
FBGN2CHROM <- file.path('../', snakemake@input[['fbgn2chrom']])
OUTDIR <- file.path('../', dirname(snakemake@output[[1]]))
REP <- snakemake@params[['rep']]
RESOLUTIONS <- snakemake@params[['resolutions']]

GENE_LOW <- snakemake@params[['gene_low']]

GENE_HIGH <- snakemake@params[['gene_high']]

UMI_LOW <- snakemake@params[['umi_low']]

UMI_HIGH <- snakemake@params[['umi_high']]

REFERENCES_DIR <- Sys.getenv('REFERENCES_DIR', '/data/LCDB/lcdb-references')



tenX.data <- Read10X(data.dir = DATA_DIR)
colnames(tenX.data) <- paste(REP, colnames(tenX.data), sep = "_")
sobj <- CreateSeuratObject(
    raw.data = tenX.data,
    min.cells = 3,
    min.genes = 200,
    project = 'GroupCluster',
    display.progress = FALSE
)
sobj@meta.data$rep <- REP
remove(tenX.data)





# Filter based on the number of genes
sobj <- FilterCells(
    object = sobj,
    subset.names = "nGene",
    low.thresholds = GENE_LOW,
    high.thresholds = GENE_HIGH
)

# Filter based on the total UMI
sobj <- FilterCells(
    object = sobj,
    subset.names = "nUMI",
    low.thresholds = UMI_LOW,
    high.thresholds = UMI_HIGH
)

### Normalization and Dimension Reduction
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

sobj <- RunPCA(
    object = sobj, 
    pc.genes = sobj@var.genes, 
    do.print = FALSE 
)

# Project other genes onto the PCA components
sobj <- ProjectPCA(object = sobj, do.print = FALSE)

### Clustering 
sobj <- FindClusters(
    object = sobj, 
    reduction.type = "pca", 
    resolution = RESOLUTIONS, 
    dims.use = 1:12, 
    save.SNN = FALSE, 
    print.output = FALSE
)
