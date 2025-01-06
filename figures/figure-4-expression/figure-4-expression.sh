#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper, Timothy R. Chu, Ali Oku

################################################################# /COPYRIGHT ###
################################################################################

## This bash script contains Figure 4 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 4a - PCA plot of gene expression after batch correction
Rscript $SRCDIR/plt-fig4a-pca-plot.r


## Figure 4b - Heatmap of differentially expressed genes 
## Figure 4c - Heatmap of most variable genes 
Rscript $SRCDIR/plt-fig4b-fig4c-heatmaps_deg_top_variable_genes.R


## Figure 4d - gTrack and junction read support heatmap of AR locus
## Left: gTrack
Rscript $SRCDIR/plt-fig4d-left-gtrack.r \
        --samples $tumor1,$tumor2,$tumor3 \
        --id_map $ID_MAP \
        --jabba_gg $gg \
        --bed $AR_BED \
        --padding 5E5 \
        --junc_table $JABBA_DIR/junction-support/WCM63-AR-footprint-union.junction_support.bed \
        --out_file $SRCDIR/fig4d-left-gtrack.pdf

## Right: Heatmap
Rscript $SRCDIR/plt-fig4d-right-junction-support-heatmap.r \
        --in_file $JABBA_DIR/junction-support/WCM63-AR-footprint-union.junction_support.bed \
        --id_map $ID_MAP \
        --metadata $METADATA \
        --out_file $SRCDIR/fig4d-right-junction-support-heatmap.pdf


## Figure 4e - Gene set enrichment results 
Rscript $SRCDIR/plt-fig4e-gsea-kegg-plots.r