#!/nfs/sw/R/R-4.0.0/bin/Rscript
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

##Load libraries 
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggfortify)


## Read library normalized count data
batch_cor_pm63 <- read.csv("Results_batch_correction/NEPC_PM63_GRCh37_counts_libnorm_batch_corrected_using_RA_cohort_minus_PM0.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

## Log transform count data
batch_cor_pm63 <- log2(batch_cor_pm63 + 1)

## Set ggplot2 themes
txt_pt2 <- 16
ttl_pt2 <- 18

txt_theme_2 <- theme(
  axis.text = element_text(size = txt_pt2),
  axis.title = element_text(size = ttl_pt2),
  plot.title = element_text(size = ttl_pt2)
)


## Compute PCA
pcaNorm <- prcomp(t(batch_cor_pm63[, 1:12]))

# Get PCA components
pcaNormComp <- data.frame(pcaNorm$x)
percent_variance_per_PC <- summary(pcaNorm)$importance[2, ] * 100

# Color PCA plot by sample type (use suffix separared by "_")
annot <- sapply(strsplit(rownames(pcaNormComp), "_"), function(x) x[[2]])

## Get histology info
pcaNormComp$histology <- sample_type_manu
sample_col <- ifelse(sample_type == "CRPC-Ad", "#52BCA3", "#CC61B0")


##########
## Plot ## 
##########

ggplot2::ggplot(pcaNormComp, ggplot2::aes(x = PC1, y = PC2, color = histology)) +
  ggplot2::geom_point(size = 5) +
  xlab(paste0("PC 1 ", "(", unname(percent_variance_per_PC[1]), "%)")) +
  ylab(paste0("PC 2 ", "(", unname(percent_variance_per_PC[2]), "%)")) +
  ggplot2::labs(color = "Histology") +
  ggplot2::scale_color_manual(values = c("#52BCA3", "#CC61B0")) +
  theme_bw() +
  txt_theme_2 +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(), 
    axis.line.x = ggplot2::element_line(color = "black"), 
    axis.line.y = ggplot2::element_line(color = "black"), 
    axis.ticks.length = grid::unit(0.3, "lines"),
    axis.ticks.x = ggplot2::element_line(color = "black"),
    axis.ticks.y = ggplot2::element_line(color = "black"),
    axis.line = ggplot2::element_line(size = 0.5),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  # Add text labels to points with ggrepel to avoid overlapping
  ggrepel::geom_text_repel(
    data = subset(pcaNormComp),
    aes(label = paste0(as.character(annot))),
    color = "black",
    hjust = 1,
    size = 4,
    segment.size = 0.0,
    nudge_x = 0.5,
    nudge_y = 0.5,
    min.segment.length = unit(0, "lines"),
    direction = "both",
    segment.alpha = 0.5,
    show.legend = FALSE
  )
