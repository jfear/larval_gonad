library(Seurat)
library(dplyr)

#' Function to read in 10x data
load_and_filter_10x <- function(data_dir, good_cells, fbgn2symbol, rep){

    # Load 10x data
    tenX.data <- Read10X(data.dir = data_dir)

    # Keep good cells
    colnames(tenX.data) <- as.vector(sapply(colnames(tenX.data), function(x) paste(rep, x, sep='_')))
    tenX.data <- tenX.data[, good_cells]

    # Rename to gene symbol
    row.names(tenX.data) <- left_join(tibble(FBgn = row.names(tenX.data)), fbgn2symbol, by='FBgn') %>% pull(gene_symbol)

    return(tenX.data)
}
