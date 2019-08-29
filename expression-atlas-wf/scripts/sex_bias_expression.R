library(DESeq2)
library(tibble)
library(dplyr)
library(ggplot2)
set.seed(42)

SAMPLE_TABLE <- snakemake@input[["sample_table"]]
COUNTS_TABLE <- snakemake@input[["counts_table"]]
OUTPUT_FILE <- snakemake@output[[1]]
SPECIES <- snakemake@wildcards[["species"]]
TISSUE <- snakemake@wildcards[["tissue"]]
LOG <- file(snakemake@log[[1]], open = "wt")
sink(LOG)
sink(LOG, type = "message")

# Debug Settings
# SAMPLE_TABLE <- "../config/sampletable.tsv"
# COUNTS_TABLE <- "../../output/expression-atlas-wf/raw_counts.feather"
# SPECIES <- "w1118"
# TISSUE <- "GO"

# Load metadata
sampletable <- read.csv(SAMPLE_TABLE, stringsAsFactors = FALSE, sep = "\t") %>%
  tidyr::separate(samplename, c("species", "tissue", "sex", "rep"), remove = FALSE) %>%
  mutate(sex = relevel(as.factor(sex), "f")) %>%
  filter(species == SPECIES, tissue == TISSUE) %>%
  column_to_rownames("samplename")

# Load count data
raw_counts <- feather::read_feather(
  COUNTS_TABLE,
  columns = c("FBgn", row.names(sampletable))
) %>%
  replace(is.na(.), 0) %>%
  filter(rowSums(select_if(., is.numeric)) != 0) %>%
  arrange(FBgn) %>%
  column_to_rownames("FBgn")

# Build deseq model
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = sampletable,
  design = ~sex
)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, betaPrior = FALSE)

# get results
res <- results(dds, name = "sex_m_vs_f", alpha = 0.01)
print(summary(res))

res.lfc <- lfcShrink(
  dds,
  coef = 2,
  res = res,
  type = "apeglm",
)

deg_res <- as.data.frame(res.lfc) %>% rownames_to_column("FBgn")

# Save results
write.table(
  deg_res,
  file = OUTPUT_FILE,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
