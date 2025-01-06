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

## This bash script contains Figure 2 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 2a - Oncoprint summarizing common prostate cancer alterations (see attached html)

## Figure 2b - Barplot of TMB and mutational signature contributions (see attached html)

## Figure 2c - Clone tree with mutational signatures superimposed on nodes (see attached html)

## Figure 2d - Barplot of tree-constrained clone proportions 
Rscript $SRCDIR/plt-fig2d-clone-barplot.r \
    --in_file $PHYCLONE_CCF_TABLE \
    --id_map $ID_MAP \
    --out_file $SRCDIR/fig2d-clone-barplot.svg