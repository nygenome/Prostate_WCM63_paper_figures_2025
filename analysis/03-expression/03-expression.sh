#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Heather Geiger

################################################################# /COPYRIGHT ###
################################################################################

## This bash script contains code associated with the gene expression analysis
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))

## Run batch correction, and DESeq2 library-size-normalization after that
Rscript $SRCDIR/batch_correction.r

#Run differential expression
Rscript $SRCDIR/differential_expression.r
