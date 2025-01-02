#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper & Timothy R. Chu

################################################################# /COPYRIGHT ###
################################################################################

## This bash script contains Figure 5 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 5a - PSA levels, treatment history, and detection rate of selected clones 

## Figure 5c - detection rate of clones shared by >1 sample 
Rscript $SRCDIR/plt-fig5c-tf-time-series.r \
        --in_file_tf $WORKING_DIR/PM63-tumor-fraction-estimates.txt \
        --in_file_tree $PHYCLONE_CCF_TABLE \
        --id_map $ID_MAP \
        --out_file $SRCDIR/fig5c-tf-time-series.svg

## Figure 5d - cosine similarity between read depth profiles 
Rscript $SRCDIR/plt-fig5d-read-depth-similarity-heatmap.r \
    --in_file $WORKING_DIR/tumor-plasma-read-depth-comparison.txt \
    --id_map $ID_MAP \
    --out_file $SRCDIR/fig5d-read-depth-similarity-heatmap.svg
