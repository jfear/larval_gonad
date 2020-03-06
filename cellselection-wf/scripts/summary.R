library(dplyr)
library(feather)
library(UpSetR)
library(grid)

RAW <- snakemake@input[["raw"]]     # list of raw barcode.tsv files
FILTERED <- snakemake@input[["filtered"]]       # list of filtered barcode.tsv files
DROPUTILS <- snakemake@input[["droputils"]]      # droputils feather file

CALLS <- snakemake@output[["cell_calls"]]       # Output a table of cell calls
SVG <- snakemake@output[["upset_plot"]]       # Output an upset plot

NAMES <- snakemake@params[[1]]           # Names of different cell ranger runs
SAMPLE <- snakemake@wildcards[["sample"]]       # Name of the sample

# Helper functions

#' Read in the 10x raw and filter cell barcodes and create a flag for which
#' GEMs are cells.
get_10x <- function(raw, filtered, name) {
    rep <- paste0("rep", gsub(".*(\\d)", "\\1", SAMPLE))

    # read in the raw data to get call GEMs
    df <- read.table(raw, stringsAsFactors = F, col.names = "barcodes")
    df$cell_id <- gsub("(\\w+)-1", paste(rep, "\\1", sep = "_"), df$barcodes)
    df[name] <- 0

    # Flag cells in filtered data
    barcodes <- read.table(filtered, stringsAsFactors = F)[, 1]
    df[df$barcodes %in% barcodes, name] <- 1
    return(df %>% select("cell_id", name))
}

#' read in the dropletutils data set and create a data.frame flagging which
#' GEMs are cells.
get_drop <- function(fname) {
    df <- read_feather(fname) %>%
        select(cell_id, droputils = is_cell) %>%
        mutate(droputils = as.integer(droputils))
    return(df)
}

# Read in 10X cell calls
tenx <- plyr::join_all(
    purrr::pmap(list(raw = RAW, filtered = FILTERED, name = NAMES), get_10x),
    by = "cell_id"
)

# Read in dropUtils cell calls
drops <- get_drop(DROPUTILS)

# Merge calls and remove always empty GEMs
df <- full_join(tenx, drops, by = c("cell_id")) %>%
    filter(rowSums(.[2:ncol(.)]) > 0) %>%
    mutate(is_cell = `cellranger3-wf` + droputils == 2)

write_feather(df, CALLS)

svg(SVG)
upset(
    df, sets = c("cellranger3-wf", "droputils"), 
    order.by = "degree", mainbar.y.label = "Cell Detected"
)
grid.text(SAMPLE, x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
dev.off()