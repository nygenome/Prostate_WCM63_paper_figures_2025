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

## Figure 3b - Violin plots of junction burden by histology and site 

## Figure 3c - Upset plots of histology-private junctions

## Figure 3d - gTrack and junction read support heatmap for a shared event
## Left: gTrack

## Right: Heatmap

## Figure 3e - Heatmap of junction burden by JaBbA event classification