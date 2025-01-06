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

## Load Libraries
library(dplyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggfortify)
library(sva)
library(pheatmap)
library(seriation)
library(dendextend)

# Function to save a heatmap as a PDF
# Arguments:
#   x: The heatmap object (e.g., from pheatmap package)
#   filename: The name of the PDF file to save the heatmap
#   width: The width of the PDF file (default is 13)
#   height: The height of the PDF file (default is 15)
save_pheatmap_pdf <- function(x, filename, width = 13, height = 15) {
  # Open a PDF device for plotting with specified filename, width, and height
  pdf(filename, width = width, height = height)

  # Set the margins for the plot
  par(mar = c(width, height, 0.5, 0.5))

  # Start a new page in the PDF
  grid::grid.newpage()

  # Create a viewport layout with specified widths and heights, adjusting for margins
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      1, 1,
      widths = grid::unit(1, "npc") - grid::unit(c(2, 2), "lines"),
      heights = grid::unit(1, "npc") - grid::unit(c(2, 2), "lines")
    )
  ))

  # Push a new viewport for the layout position
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

  # Draw the heatmap gtable object within the layout
  grid::grid.draw(x$gtable)

  # Close the PDF device
  dev.off()
}


# Read lib normalized count data
batch_cor_pm63 <- read.csv("/gpfs/commons/projects/nepc/analysis/RNA/Results_batch_correction/NEPC_PM63_GRCh37_counts_libnorm_batch_corrected_using_RA_cohort_minus_PM0.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# log transform
batch_cor_pm63 <- log2(batch_cor_pm63 + 1)


# Create column corresponding to standard deviation
batch_cor_pm63$std <- data.frame(apply(batch_cor_pm63, 1, sd, na.rm = TRUE))[, 1]
batch_cor_pm63_ord <- batch_cor_pm63[order(batch_cor_pm63$std, decreasing = TRUE), ]

# Load clinical information including patient ID and manuscript ID
manu_id <- read.table(
  file = "/gpfs/commons/projects/nepc/analysis/metadata/PM63_clinical_data_sync_NYGC.tsv", 
  sep = '\t', 
  header = TRUE,
  check.names = FALSE
)

# Filter and clean the data
manu_id <- manu_id %>%
  select(`Associated sample`, `New sample name`) %>%
  filter(`Associated sample` != "") %>%
  na.omit()

# Create manuID column with modified values
manu_id <- manu_id %>%
  mutate(
    manuID = substr(`Associated sample`, 1, 8),
    manuID = sub("_$", "", manuID),
    manuID = ifelse(manuID == "PM63_Z19", "PM63_Z18", manuID)
  )

# Create a named vector mapping patient ID to manuscript ID
sample_manu_map <- setNames(manu_id$`New sample name`, manu_id$manuID)

# Extract the relevant manuIDs from batch_cor_pm63 column names
manu_ids <- unname(sample_manu_map[colnames(batch_cor_pm63)[1:12]])


# Create a named vector mapping sample IDs to their corresponding cancer type (Squamous or Adeno)
sample_type_map <- setNames(
  c("Squamous", "Squamous", "Squamous", "Squamous", "Squamous", "Squamous", 
    "Adeno", "Adeno", "Adeno", "Adeno", "Adeno", "Adeno"),
  c("PM63_Z1", "PM63_Z6", "PM63_Z7", "PM63_Z8", "PM63_Z9", "PM63_Z11", 
    "PM63_Z2", "PM63_Z3", "PM63_Z4", "PM63_Z5", "PM63_Z10", "PM63_Z18")
)

# Extract the sample type based on the column names of batch_cor_pm63
sample_type <- unname(sample_type_map[colnames(batch_cor_pm63)])

# Map sample types to a more descriptive label for adenocarcinoma and squamous cell carcinoma
sample_type_manu <- ifelse(sample_type == "Adeno", "CRPC-Ad", "CRPC-SCC")

# Convert the sample type to a factor with specified levels for consistent ordering
sample_type <- factor(sample_type_manu, levels = c("CRPC-Ad", "CRPC-SCC"))

#Use manuscript IDS
colnames(batch_cor_pm63) <- manu_ids


#-------------------------------------------------------------------------------
#----------------------- Top Variable Genes  -----------------------------------


# Prep top annotation
annotation_col <- data.frame(
  Histology = sample_type_manu,
  row.names = colnames(batch_cor_pm63[1:12])
)

# Scale heatmap colors
default_pheatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

heatmap_col <- c(
  default_pheatmap_colors[1],
  default_pheatmap_colors,
  default_pheatmap_colors[100]
)

breaksList <- c(-6, seq(from = -2, to = 2, length.out = 100), 6)

phtmap <- pheatmap::pheatmap(batch_cor_pm63_ord[0:4889, 0:12],
  show_rownames = FALSE,
  clustering_method = "complete",
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = annotation_col,
  color = heatmap_col,
  breaks = breaksList,
  main = "Top 4889 Most Variable Genes (by stdev) using batch-corrected \n Log2(libnorm, batch-cor), scaled"
)

col_dend <- phtmap[[2]]
col_dend <- dendextend::rotate(col_dend,
  order = rev(names(batch_cor_pm63_ord[0:4889, 0:12])[seriation::get_order(col_dend)])
)


# Plot heatmap of top 4889 variable genes
ann_colors <- list(
  Histology = c(`CRPC-Ad` = "#52BCA3", `CRPC-SCC` = "#CC61B0")
)

top_variable_heatmap <- pheatmap::pheatmap(batch_cor_pm63_ord[0:4889, 0:12],
  show_rownames = FALSE,
  clustering_method = "complete",
  cluster_cols = as.hclust(col_dend),
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_col,
  breaks = breaksList,
  main = "Top 4889 Most Variable Genes using batch-corrected normalized \n gene expression data (row scaled)"
)

# Save heatmap
save_pheatmap_pdf(
  x = top_variable_heatmap,
  filename = "plots/top_variable_genes_heatmap.pdf",
  width = 8, height = 8
)

#------------------------------------------------------------------------------
#------------ Differential Gene Expression (DEG) ------------------------------

# Load deg results (analysis done by Heather Geiger)
deg_results <- read.csv("DE_PM63_squam_vs_adeno_after_batch_correct.csv",
  check.names = FALSE, row.names = 1
)

deg_results <- na.omit(deg_results)

# Limit to only genes with p adj < 0.01
deg_results_sig <- deg_results[deg_results$padj < 0.01, ]

# Prep top annotation
annotation_col <- data.frame(
  Histology = sample_type_manu,
  row.names = colnames(batch_cor_pm63)[1:12]
)

#Heatmap (DEG: Squamous vs Adeno)
sig_dge_heatmap <- pheatmap::pheatmap(as.matrix(batch_cor_pm63[rownames(deg_results_sig), 1:12]),
                                      show_rownames = FALSE,
                                      clustering_method = "complete",
                                      scale = "row",
                                      annotation_col = annotation_col,
                                      annotation_colors = ann_colors,
                                      color = heatmap_col,
                                      breaks = breaksList,
                                      main = "DE squam vs. adeno using batch corrected, n=4889, FDR 1% \n Log2(libnorm, batch-cor), scaled"
)

# Save heatmap (DEG: Squam vs Adeno)
save_pheatmap_pdf(
  x = sig_dge_heatmap,
  filename = "plots/sig_deg_heatmap.pdf",
  width = 8, height = 8
)


