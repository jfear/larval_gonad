library(DropletUtils)
library(dplyr)
library(feather)
library(ggplot2)
theme_set(theme_bw())

DATA_DIR <- dirname(snakemake@input[[1]])
PVAL_PLOT <- snakemake@output[["pval_plot"]]
BC_RANK_PLOT <- snakemake@output[["barcode_rank"]]
CELL_CALLS <- snakemake@output[["cell_calls"]]
SAMPLE <- snakemake@wildcards[["sample"]]
REP <- gsub(".*(\\d)", "rep\\1", SAMPLE)

# Load 10x data
tenx <- read10xCounts(DATA_DIR)

# Determine Barcode Ranks
barcode_rank <- barcodeRanks(counts(tenx))

# Estimate Empty Droplets
set.seed(42)
empty_drops <- emptyDrops(counts(tenx))
is_cell <- empty_drops$FDR <= 0.01
is_cell[is.na(is_cell)] <- FALSE

# Plot pval
ggplot(
    as.data.frame(empty_drops),
    aes(x = Total, y = -LogProb, color = is_cell)
) +
    geom_point() +
    scale_color_manual(values = c("lightgray", "black")) +
    labs(x = "Total UMI", y = "-Log P-Value", title = SAMPLE) +
    scale_x_log10() +
    scale_y_log10()

ggsave(PVAL_PLOT)

# Output Cell Selection Calls
df <- as.data.frame(barcode_rank)
df$is_cell <- is_cell
df$barcode <- tenx@colData$Barcode
df$cell_id <- gsub("(.*?)-1", paste0(REP, "_\\1"), df$barcode)
df <- df %>% select(cell_id, barcode, everything())
write_feather(df, CELL_CALLS)

# plot barcod rank
ggplot(df, aes(x = rank, y = total, color = is_cell)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    geom_hline(
        yintercept = barcode_rank$knee,
        linetype = 2, color = "dodgerBlue"
    ) +
    geom_hline(
        yintercept = barcode_rank$inflection,
        linetype = 2, color = "forestgreen"
    ) +
    scale_color_manual(values = c("lightgray", "black")) +
    labs(
        x = "Ranks", y = "Total UMI",
        title = paste0(SAMPLE, " (DropletUtils)")
    )
ggsave(BC_RANK_PLOT)