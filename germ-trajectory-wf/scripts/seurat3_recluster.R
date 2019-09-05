library(seurat3)
library(dplyr)
library(ggplot2)

if (Sys.getenv("SNAKE_DEBUG", FALSE)) {
    source("../lib/snakemake.R")
    snakemake <- snakemake_debug(
        input = list(
            robj = "../output/seurat3-cluster-wf/combined_n3.Robj",
            lit_genes = "../config/literature_genes.yaml",
            gene_annot = "../references/gene_annotation_dmel_r6-26.feather"
        ),
        params = list(germ_clusters = c(6, 4, 2, 0))
    )
}

# Load Metadata
lit_genes <- yaml::yaml.load_file(snakemake@input[["lit_genes"]])
fbgn2symbol <- feather::read_feather(snakemake@input[["gene_annot"]], columns = c("FBgn", "gene_symbol"))

# Load combined robj
load(snakemake@input[["robj"]])

# Subset out the germline clusters
sobj <- subset(combined, ident = snakemake@params[["germ_clusters"]])

# Re-process with subset
sobj <- RunPCA(sobj, npcs = 30, features = VariableFeatures(sobj), verbose = FALSE)
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:9)
sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:9, verbose = FALSE)
sobj <- FindClusters(sobj, resolution = 0.2, verbose = FALSE)

# Plot new clusters
DimPlot(sobj, reduction = "umap")

# Plot gene sets
gonia <- fbgn2symbol %>%
    filter(FBgn %in% lit_genes[["gonia"]]) %>%
    pull(gene_symbol)
FeaturePlot(sobj, features = gonia)

ps <- fbgn2symbol %>%
    filter(FBgn %in% lit_genes[["spermatocytes"]]) %>%
    pull(gene_symbol)
FeaturePlot(sobj, features = ps)

# Look at new biomarkers
biomarkers <- FindAllMarkers(sobj, only.pos = TRUE, verbose = FALSE)
biomarkers %>%
    filter(cluster == 0) %>%
    arrange(gene) %>%
    pull(gene)
