library(Seurat)

# Get snakemake resources
WORKDIR <- file.path('../', dirname(snakemake@input[[1]]))
RESOLUTIONS <- snakemake@params[['resolutions']]

REFERENCES_DIR <- Sys.getenv('REFERENCES_DIR', '/data/LCDB/lcdb-references')
fbgn2gene <- read.table(
    file.path(REFERENCES_DIR, 'dmel/r6-24/fb_annotation/dmel_r6-24.fb_annotation'),
    sep='\t',
    stringsAsFactors=F,
    header=T,
    quote=NULL,
    fill=TRUE,
    colClasses=c("character", "NULL", "character", rep("NULL", 6))
) %>% distinct()

# Load the Seurat Object
load(file.path(WORKDIR, 'seurat.Robj'))
sobj <- object

# Run differential expression
for (i in RESOLUTIONS){
  name <- paste0('res.', i)
  sobj <- SetAllIdent(sobj, id = name)

  # Get markers ignoring replicate.
  fname <- paste0('biomarkers_', name, '.tsv')
  markers <- FindAllMarkers(object = sobj, only.pos = TRUE, print.bar = FALSE)
  markers = merge(fbgn2symbol, markers, by.x='primary_FBgn', by.y='gene', all.y=T)
  save_biomarkers(markers = markers, dir = OUTDIR, fname = fname)
}
