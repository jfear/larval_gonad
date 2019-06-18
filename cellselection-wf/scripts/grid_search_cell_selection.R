library(Seurat)
library(dplyr)
library(feather)

SAMPLE <- snakemake@wildcards[['sample']]
DATA_DIR <- dirname(snakemake@input[[1]])
OUTPUT <- snakemake@output[[1]]

# Build list of combinations for Grid Search
GENE_LOW <- c(200, 500, 1000)
GENE_HIGH <- c(4000, 5000, 6000, Inf)
PERCENT_MITO <- c(20, 40, 60, Inf)
PERCENT_RIBO <- c(2, 4, 5, Inf)
COMBOS <- expand.grid(GENE_LOW, GENE_HIGH, PERCENT_MITO, PERCENT_RIBO)

# Basic Seurat Clustering using 12 PCs
run <- function(tenX.data, gene_low, gene_high, pct_mito, pct_ribo){
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
        res = paste0('res.', RESOLUTIONS[[i]])
        lvls <- metadata[res] %>% distinct() %>% pull()
        if (length(lvls) >= 13){
           df <- metadata %>% 
               select("cell_id", "nGene", "nUMI", percent_mito = "percent.mito", percent_ribo = "percent.ribo", "gene_low_threshold", "gene_high_threshold", "percent_mito_threshold", "percent_ribo_threshold", cluster = res)
           return(df)
        }
    }
}

# Run Grid Search
tenX.data <- Read10X(data.dir = DATA_DIR)
results <- list()
for (i in 1:dim(COMBOS)[[1]]){
    cbn <- COMBOS[1, ]
    results[[i]] <- run(tenX.data, cbn[[1]], cbn[[2]], cbn[[3]], cbn[[4]])
}

# Write out table with all results
df <- bind_rows(results)
write_feather(df, OUTPUT)