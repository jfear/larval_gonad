library(tibble)
library(dplyr)
library(ggplot2)
library(cowplot)

GVCO <- "/home/fearjm/Downloads/neo-X_pseudo/gonia_vs_cytes(1).tsv"
GVC <- "../../output/seurat3-cluster-wf/combined_n3_gonia_vs_cytes.feather"
MULLERE <- "../../output/neox-wf/sturgill_2007_mullerE.feather"
MULLERD <- "../../output/neox-wf/sturgill_2007_mullerD.feather"

gvc <- as_tibble(read.csv(GVCO, header = TRUE, sep = "\t", stringsAsFactor = FALSE)) %>%
    rename(FBgn = primary_FBgn, pct_cytes = `pct.2`) %>%
    select(FBgn, pct_cytes)

md <- feather::read_feather(MULLERD) %>%
    rename(mullerD = `Gene fate`, FBgn = `D.mel`) %>%
    left_join(gvc) %>%
    filter(!is.na(pct_cytes))

me <- feather::read_feather(MULLERE) %>%
    rename(mullerE = `Fate`, FBgn = `D.mel`) %>%
    left_join(gvc) %>%
    filter(!is.na(pct_cytes))

p1 <- ggplot(md, aes(x = mullerD, y = pct_cytes)) +
    geom_boxplot() + ggtitle("Muller D") + 
    theme(axis.title.x = element_blank())

p2 <- ggplot(me, aes(x = mullerE, y = pct_cytes)) +
    geom_boxplot() + ggtitle("Muller E") +
    theme(axis.title.x = element_blank())

plot_grid(p1, p2)


# Muller E
wilcox.test(
    me %>% filter(mullerE == "Conserved" | mullerE == "Move On") %>% pull(pct_cytes),
    me %>% filter(mullerE == "Gene Death" | mullerE == "Move Off") %>% pull(pct_cytes)
)

# Muller D
wilcox.test(
    md %>% filter(mullerD == "Conserved" | mullerD == "Move On") %>% pull(pct_cytes),
    md %>% filter(mullerD == "Gene Death" | mullerD == "Move Off") %>% pull(pct_cytes)
)