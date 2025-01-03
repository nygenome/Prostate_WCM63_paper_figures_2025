#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Timothy R. Chu

################################################################# /COPYRIGHT ###
################################################################################

## This bash script contains code associated with the mutational signatures analysis
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


#####################
## deconstructSigs ##
#####################
./SigMatrixGen.py -o $project -r GRCh38 -i $DECONSTRUCTSIGS_DIR
Rscript run_deconstructSigs.R -d $DECONSTRUCTSIGS_DIR --highconf -o $output_sbs -v v3.2
Rscript run_deconstructSigs_DBS.R -d $DECONSTRUCTSIGS_DIR --highconf -o $output_dbs -v v3.2
Rscript run_deconstructSigs_ID.R -i $DECONSTRUCTSIGS_DIR/output/ID/$project.ID83.all -o $output_id -v 3.2

