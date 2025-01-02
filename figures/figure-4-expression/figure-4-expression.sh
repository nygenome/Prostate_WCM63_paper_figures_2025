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

## Figure 4b - Heatmap of differentially expressed genes 

## Figure 4c - Heatmap of most variable genes 

## Figure 4d - gTrack and junction read support heatmap of AR locus
## Left: gTrack

## Right: Heatmap

## Figure 4e - Gene set enrichment results 
