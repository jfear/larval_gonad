library(Seurat)
library(dplyr)
source('../lib/seurat.R')

SAMPLE <- snakemake@wildcards[['sample']]
WORKDIR <- dirname(snakemake@input[[1]])
load(file.path(WORKDIR, 'seurat.Robj'))

if (SAMPLE == 'combined'){
    sobj <- combined
}

dump_seurat(sobj, WORKDIR)
