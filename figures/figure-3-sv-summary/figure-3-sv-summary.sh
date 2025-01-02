#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper 

################################################################# /COPYRIGHT ###
################################################################################

## This bash script contains Figure 3 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 3a - Heatmap showing JaBbA event footprints across the genome 
Rscript $SRCDIR/plt-fig3a-event-footprint-heatmap.r \
    --in_dir $JABBA_DIR/jabba \
    --junction_support_mtx $JABBA_DIR/junction-support/PM63.all_junctions.junction_support.bed \
    --tn_file $TNFILE \
    --id_map $ID_MAP \
    --chr_len $CHR_LEN \
    --metadata $METADATA \
    --binsize 1E6 \
    --out_file $SRCDIR/fig3a-event-footprint-heatmap.pdf


## Figure 3b - Violin plots of junction burden by histology and site 
Rscript $SRCDIR/plt-fig3b-junction-burden-violin.r \
    --in_file $JABBA_DIR/jabba.events-junction-summary.csv \
    --tn_file $TNFILE \
    --metadata $METADATA \
    --out_file $SRCDIR/fig3b-junction-burden-violin.svg


## Figure 3c - Upset plots of histology-private junctions
Rscript $SRCDIR/plt-fig3c-upset-histology-private-junctions.r \
        --in_file $JABBA_DIR/junction-support/PM63.all_junctions.junction_support.bed \
        --id_map $ID_MAP \
        --tn_file $TNFILE \
        --metadata $METADATA \
        --out_file $SRCDIR/fig3c-upset-histology-private-junctions.pdf


## Figure 3d - gTrack and junction read support heatmap for a shared event
## Left: gTrack
Rscript $SRCDIR/plt-fig3d-left-gtrack.r \
        --samples $tumor1,$tumor2,$tumor3 \
        --id_map $ID_MAP \
        --jabba_gg $gg \
        --bed $JABBA_DIR/junction-support/PM63-chr5-chromothripsis-footprint-union.bed \
        --padding 5E5 \
        --junc_table $JABBA_DIR/junction-support/PM63-chr5-chromothripsis-footprint-union.junction_support.bed \
        --out_file $SRCDIR/fig3d-left-gtrack.pdf

## Right: Heatmap
Rscript $SRCDIR/plt-fig3d-right-junction-support-heatmap.r \
        --in_file $JABBA_DIR/junction-support/PM63-chr5-chromothripsis-footprint-union.junction_support.bed \
        --id_map $ID_MAP \
        --metadata $METADATA \
        --out_file $SRCDIR/fig3d-right-junction-support-heatmap.pdf


## Figure 3e - Heatmap of junction burden by JaBbA event classification
Rscript $SRCDIR/plt-fig3e-junction-burden-heatmap.r \
    -k 2 \
    --in_file $JABBA_DIR/jabba.events-junction-summary.csv \
    --tn_file $TNFILE \
    --metadata $METADATA \
    --id_map $ID_MAP \
    --out_file $SRCDIR/fig3e-junction-burden-heatmap.pdf
    