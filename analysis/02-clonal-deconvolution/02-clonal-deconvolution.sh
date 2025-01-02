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

## This bash script contains code associated with the clonal deconvolution / phylogenetic analysis
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Generate conipher input (only ever used for PyClone-VI)


################
## PyClone-VI ##
################

Rscript $SRCDIR/init-convert-conipher-input.r \
        -i $ANALYSIS_DIR/conipher/conipher_input_fixed.txt \
        -g male \
        -o $WORKING_DIR/WCM63.pyclone_vi_input.conipher_converted.txt

pyclone-vi fit \
        -i $WORKING_DIR/WCM63.pyclone_vi_input.conipher_converted.txt \
        -o $WORKING_DIR/WCM63.pyclone_vi_output.100.20240906.h5 \
        -c 40 \
        -d beta-binomial \
        -r 100 \
        --num-threads 4 \
        --seed 20240906 

pyclone-vi write-results-file \
                -i $WORKING_DIR/PM63.pyclone_vi_output.100.20240906.h5 \
                -o $WORKING_DIR/PM63.pyclone_vi_output.100.20240906.txt



##############
## Phyclone ##
##############
phyclone run \
    -i $in_file \
    -c $in_file_clusters \
    --burnin 3000 \
    --num-iters 50000 \
    --num-chains 16 \
    --seed 20240920 \
    -o $out_file_pkl

phyclone map \
    -i $out_file_pkl \
    -t $out_file_nwk \
    -o $out_file_tsv

Rscript $SRCDIR/init-post-process-phyclone.r \
    --in_file_tsv $out_file_tsv \
    --in_file_nwk $out_file_nwk \
    --out_file_pdf $out_file_tree_pdf \
    --out_file_txt $out_file_txt
