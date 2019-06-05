library(Seurat)
library(dplyr)
source('../lib/seurat.R')

# Get snakemake resources
WORKDIR <- dirname(snakemake@input[[1]])
RESOLUTION <- snakemake@wildcards[['resolution']]

REFERENCES_DIR <- Sys.getenv('REFERENCES_DIR', '/data/LCDB/lcdb-references')
fbgn2symbol <- read.table(
    file.path(REFERENCES_DIR, 'dmel/r6-24/fb_annotation/dmel_r6-24.fb_annotation'),
    sep = '\t',
    stringsAsFactors = FALSE,
    header = TRUE,
    quote = NULL,
    fill = TRUE,
    colClasses = c("character", "NULL", "character", rep("NULL", 6))
) %>% distinct()

# Load the Seurat Object
load(file.path(WORKDIR, 'seurat.Robj'))

if (snakemake@wildcards[['sample']] == 'combined'){
  sobj <- combined
}

# Run differential expression
name <- paste0('res.', RESOLUTION)
sobj <- SetAllIdent(sobj, id = name)

# Get markers 
fname <- paste0('biomarkers_', name, '.tsv')
markers <- FindAllMarkers(object = sobj, only.pos = TRUE, print.bar = FALSE)
markers <- merge(fbgn2symbol, markers, by.x = 'primary_FBgn', by.y = 'gene', all.y = TRUE)
save_biomarkers(markers = markers, dir = WORKDIR, fname = fname)
