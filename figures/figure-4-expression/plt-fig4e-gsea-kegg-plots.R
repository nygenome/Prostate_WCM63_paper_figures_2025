################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Ali E. Oku

################################################################# /COPYRIGHT ###
################################################################################
## Visualize differential gene expression using heatmaps


## Load libraries
library(dplyr)
library(ggplot2)



## Squamous vs Adeno (GO - BP)
## Read GSEA results from webgestalt
gsea_data_mf <- read.table("data/de_squamous_vs_adeno_gsea_go_mf.txt", 
                           sep = "\t", 
                           header = TRUE, check.names = FALSE)

# Make description all upper case
gsea_data_mf$description <- toupper(gsea_data_mf$description)

###############
## Plot GSEA ##
###############
ggplot(gsea_data_mf, aes(normalizedEnrichmentScore, fct_reorder(description, normalizedEnrichmentScore)), showCategory=(n*2)) + 
  ggstance::geom_barh(stat='identity', fill="#CC61B0") + 
  theme_bw() + txt_theme_2 + ylab(NULL) + xlab("NES") +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(color = "black"),
    axis.line.y = ggplot2::element_line(color = "black"),
    axis.ticks.length = grid::unit(0.3, "lines"),
    axis.ticks.x = ggplot2::element_line(color = "black"),
    axis.ticks.y = ggplot2::element_line(color = "black"),
    axis.line = ggplot2::element_line(size = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(filename = "plots/gsea_plot_squamous_vs_adeno_go_mf.pdf", 
                width = 12, 
                height = 8, 
                dpi = 320)


## Squamous vs Adeno (KEGG)
## Read KEGG results from webgestalt
gsea_data_kegg <- read.table("data/de_squamous_vs_adeno_gsea_kegg.txt", 
                             sep = "\t", 
                             header = TRUE, 
                             check.names = FALSE)

gsea_data_kegg$description <- toupper(gsea_data_kegg$description)
gsea_data_kegg$nesColor <- ifelse(gsea_data_kegg$normalizedEnrichmentScore > 0, "#CC61B0", "#52BCA3")

###############
## Plot KEGG ##
###############
ggplot(gsea_data_kegg, aes(normalizedEnrichmentScore, 
                           fct_reorder(description, normalizedEnrichmentScore)), showCategory=(n*2)) + 
  ggstance::geom_barh(stat='identity', fill=gsea_data_kegg$nesColor) + 
  theme_bw() + txt_theme_2 + ylab(NULL) + xlab("NES") +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),  # Remove all borders
    axis.line.x = ggplot2::element_line(color = "black"),
    axis.line.y = ggplot2::element_line(color = "black"),
    axis.ticks.length = grid::unit(0.3, "lines"),
    axis.ticks.x = ggplot2::element_line(color = "black"),
    axis.ticks.y = ggplot2::element_line(color = "black"),
    axis.line = ggplot2::element_line(size = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(filename = "plots/gsea_plot_squamous_vs_adeno_kegg.pdf", 
                width = 14, 
                height = 8, 
                dpi = 320)

