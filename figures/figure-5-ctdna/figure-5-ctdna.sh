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

## This bash script contains Figure 5 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 5a - PSA levels, treatment history, and detection rate of selected clones 

## Figure 5c - detection rate of clones shared by >1 sample 

## Figure 5d - cosine similarity between read depth profiles 