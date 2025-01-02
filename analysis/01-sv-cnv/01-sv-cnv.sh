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

## This bash script contains code associated with the structural variation analysis
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))



##########################################
## Generate JaBbA profiles and SV calls ##
########################################## 

## Run fragcounter to get read depth profiles in 1KB bins 
bash $SRCDIR/jabba/preprocessing.sh fragcounter $SRCDIR/jabba_config.txt

## Generate dryclean PON from normal read depth profiles
bash $SRCDIR/jabba/preprocessing.sh dryclean-pon $SRCDIR/jabba_config.txt

## De-noise tumor read depth profiles with dryclean 
bash $SRCDIR/jabba/preprocessing.sh dryclean $SRCDIR/jabba_config.txt

## Run JaBbA to infer junction-balanced genome graphs
bash $SRCDIR/jabba/jabba.sh $SRCDIR/jabba_config.txt

## Run calling to infer simple/complex SVs
bash $SRCDIR/jabba/calling.sh $SRCDIR/jabba_config.txt



###########################################################
## Junction read support of selected events for plotting ##
###########################################################

## chr5 chromothripsis 
Rscript $SRCDIR/init-footprint-union.r \
    --in_dir $JABBA_DIR/jabba \
    --tn_file $TNFILE \
    --chr chr5 \
    --event_type chromothripsis \
    --out_file $JABBA_DIR/junction-support/WCM63-chr5-chromothripsis-footprint-union.bed

Rscript $SRCDIR/patient-junction-support-by-bed.r \
    --jba_dir $JABBA_DIR/jabba \
    --bam_dir $PROJECT_DIR/data \
    --tn_file $TNFILE \
    --bed $JABBA_DIR/junction-support/WCM63-chr5-chromothripsis-footprint-union.bed \
    --genome $GENOME \
    --out_file $JABBA_DIR/junction-support/WCM63-chr5-chromothripsis-footprint-union.junction_support.bed


## chr2 chromoplexy 
Rscript $SRCDIR/init-footprint-union.r \
    --in_dir $JABBA_DIR/jabba \
    --tn_file $TNFILE \
    --chr chr2 \
    --event_type chromoplexy \
    --out_file $JABBA_DIR/junction-support/WCM63-chr2-chromoplexy-footprint-union.bed

Rscript $SRCDIR/patient-junction-support-by-bed.r \
    --jba_dir $JABBA_DIR/jabba \
    --bam_dir $PROJECT_DIR/data \
    --tn_file $TNFILE \
    --bed $JABBA_DIR/junction-support/WCM63-chr2-chromoplexy-footprint-union.bed \
    --genome $GENOME \
    --out_file $JABBA_DIR/junction-support/WCM63-chr2-chromoplexy-footprint-union.junction_support.bed


## chr3 bfb 
Rscript $SRCDIR/init-footprint-union.r \
    --in_dir $JABBA_DIR/jabba \
    --tn_file $TNFILE \
    --chr chr3 \
    --event_type bfb \
    --out_file $JABBA_DIR/junction-support/WCM63-chr3-bfb-footprint-union.bed

Rscript $SRCDIR/patient-junction-support-by-bed.r \
    --jba_dir $JABBA_DIR/jabba \
    --bam_dir $PROJECT_DIR/data \
    --tn_file $TNFILE \
    --bed $JABBA_DIR/junction-support/WCM63-chr3-bfb-footprint-union.bed \
    --genome $GENOME \
    --out_file $JABBA_DIR/junction-support/WCM63-chr3-bfb-footprint-union.junction_support.bed


## AR 
Rscript $SRCDIR/patient-junction-support-by-bed.r \
        --jba_dir $JABBA_DIR/jabba \
        --bam_dir $PROJECT_DIR/data \
        --tn_file $TNFILE \
        --bed $AR_BED \
        --pad 1E5 \
        --genome $GENOME \
        --out_file $JABBA_DIR/junction-support/WCM63-AR-footprint-union.junction_support.bed



#####################################
## Junction support across all SVs ##
##################################### 

Rscript $SRCDIR/patient-junction-support-by-bed.r \
        --jba_dir $JABBA_DIR/jabba \
        --bam_dir $PROJECT_DIR/data \
        --tn_file $TNFILE \
        --genome $GENOME \
        --out_file $JABBA_DIR/junction-support/WCM63.all_junctions.junction_support.bed


## Summarize junction conservation within events 
Rscript $SRCDIR/plt-event-union-junction-support-heatmap.r \
        --in_file $JABBA_DIR/junction-support/all-events/WCM63-event-union-junction-support.txt \
        --metadata $METADATA \
        --id_map $ID_MAP \
        --out_file $JABBA_DIR/junction-support/all-events/fig/WCM63-event-union-junction-support.pdf


## Compare alterations to DE genes 
Rscript $SRCDIR/init-compare-alterations-to-de-genes.r \
        --onco_mtx $ANALYSIS_DIR/oncoprint/oncomatrix_WCM63.txt \
        --de_genes $DE_GENES \
        --metadata $METADATA
