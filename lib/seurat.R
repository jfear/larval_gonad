# Set of useful R function for use with seurat.
#' Saves data from the seurat object for use by other tools.
#'
#' @param object Seurat object
#' @param dir Output directory
dump_seurat <- function(object, dir) {
    meta_data <- as.data.frame(object@meta.data)
    write.table(meta_data, file = file.path(dir, 'metadata.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

    raw_data <- as.data.frame(as.matrix(object@raw.data))
    write.table(raw_data, file = file.path(dir, 'raw.tsv'), quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)

    norm <- as.data.frame(as.matrix(object@data))
    write.table(norm, file = file.path(dir, 'normalized_read_counts.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

    if (length(object@hvg.info) > 0){
	    dispersion <- as.data.frame(as.matrix(object@hvg.info))
	    write.table(dispersion, file = file.path(dir, 'dispersion.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    }

    var_genes = object@var.genes
    write(var_genes, file = file.path(dir, 'var_genes.txt'))

    scaled <- as.data.frame(object@scale.data)
    write.table(scaled, file = file.path(dir, 'scaled.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

    if (length(object@dr$pca) > 0){
	    pca_res <- as.data.frame(object@dr$pca@cell.embeddings)
	    write.table(pca_res, file = file.path(dir, 'principal_components_cell.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

	    gene_loadings <- as.data.frame(object@dr$pca@gene.loadings)
	    write.table(gene_loadings, file = file.path(dir, 'principal_components_gene.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

	    pca_stdev <- object@dr$pca@sdev
	    write.table(pca_stdev, file = file.path(dir, 'principal_components_stdev.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    }

    if (length(object@dr$cca) > 0){
	    pca_res <- as.data.frame(object@dr$cca@cell.embeddings)
	    write.table(pca_res, file = file.path(dir, 'cca_cell.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

	    gene_loadings <- as.data.frame(object@dr$cca@gene.loadings)
	    write.table(gene_loadings, file = file.path(dir, 'cca_gene.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

	    gene_loadings_full <- as.data.frame(object@dr$cca@gene.loadings.full)
	    write.table(gene_loadings_full, file = file.path(dir, 'cca_gene_full.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    }

    if (length(object@dr$cca.aligned) > 0){
	    pca_res <- as.data.frame(object@dr$cca.aligned@cell.embeddings)
	    write.table(pca_res, file = file.path(dir, 'cca_aligned_cell.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    }

    clusters <- object@meta.data[, grepl('res', colnames(object@meta.data))]
    write.table(clusters, file = file.path(dir, 'clusters.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

    tsne <- as.data.frame(object@dr$tsne@cell.embeddings)
    write.table(tsne, file = file.path(dir, 'tsne.tsv'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

}

#' Saves biomarker table.
#'
#' @param markers biomarkers table
#' @param dir Output directory
#' @param fname Names to output
save_biomarkers <- function(markers, dir, fname = 'biomarkers.tsv') {
	write.table(markers, file = file.path(dir, fname), quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}
