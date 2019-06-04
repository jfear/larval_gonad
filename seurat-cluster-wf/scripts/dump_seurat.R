library(dplyr)
source('../../lib/seurat.R')

FNAME <- snakemake@input[[1]]
OUTDIR <- dirname(FNAME)

load(FNAME)

dump_seurat(object, OUTDIR)