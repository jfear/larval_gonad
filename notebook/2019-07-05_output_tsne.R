library(Seurat)
library(dplyr)
setwd("notebook")

load("../output/seurat3-cluster-wf/combined_n3.Robj")

combined <- RunTSNE(combined, reduction = "pca", dims = 1:9)
tsne <- as_tibble(as.data.frame(Embeddings(combined, reduction = "tsne")) %>%
    tibble::rownames_to_column("cell_id"))

feather::write_feather(tsne, "../output/notebook/2019-07-05_output_tsne.feather")
