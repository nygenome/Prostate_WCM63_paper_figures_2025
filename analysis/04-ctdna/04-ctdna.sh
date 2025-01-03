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

## This bash script contains code associated with the ctDNA analysis
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))



###############################################################
## Run pileup on plasma samples with union of tumor variants ##
###############################################################

## WCM63
while read tumor normal; do

    tn=$tumor--$normal
    tumor_bam=$CFDNA_DIR/Sample_$tumor/compbio/analysis/$tumor.final.bam
    normal_bam=$CFDNA_DIR/Sample_$normal/analysis/$normal.final.bam

    dummy_vcf=$WORKING_DIR/pileup/$tn.dummy.vcf
    out_vcf=$WORKING_DIR/pileup/$tn.union.vcf

    ## We need to check all variants, so create a dummy vcf with no variants
    vcf_header="##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$normal\t$tumor"
    grep '^#' $WCM63_UNION_VCF | head -n -1 > $dummy_vcf
    echo -e $vcf_header >> $dummy_vcf

    bgzip -f $dummy_vcf
    tabix $dummy_vcf.gz

    python $SRCDIR/get_vaf.py \
        --tumor_bam $tumor_bam \
        --normal_bam $normal_bam \
        --tumor $tumor \
        --normal $normal \
        --vcf $dummy_vcf.gz \
        --union-vcf $WCM63_UNION_VCF \
        --output $out_vcf

done < $TNFILE


## Controls
while read tumor normal; do

    tn=$tumor--$normal

    tumor_bam=$CFDNA_DIR/Sample_$tumor/compbio/analysis/$tumor.final.bam
    normal_bam=$CFDNA_DIR/Sample_$normal/analysis/$normal.final.bam

    dummy_vcf=$WORKING_DIR/control/$tn.dummy.vcf
    out_vcf=$WORKING_DIR/control/$tn.union.vcf

    # ## We need to check all variants, so create a dummy vcf with no variants
    vcf_header="##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$normal\t$tumor"
    grep '^#' $WCM63_UNION_VCF | head -n -1 > $dummy_vcf
    echo -e $vcf_header >> $dummy_vcf

    bgzip -f $dummy_vcf
    tabix $dummy_vcf.gz

    python $SRCDIR/get_vaf.py --tumor_bam $tumor_bam \
        --normal_bam $normal_bam \
        --tumor $tumor \
        --normal $normal \
        --vcf $dummy_vcf.gz \
        --union-vcf $WCM63_UNION_VCF \
        --output $out_vcf

done < $CONTROLS



##############################################################
## Compute per-clone detection rate, overall detection rate ##
##############################################################

## WCM63
while read tumor normal; do
    
        tn=$tumor--$normal
        plasma_vcf=$WORKING_DIR/pileup/$tn.union.vcf
        out_file=$WORKING_DIR/summary/$tumor.detection_summary.txt

        Rscript $SRCDIR/init-detection-summary.r \
                --clones $CLONE_ASSIGNMENTS \
                --vcf $plasma_vcf \
                --sample_id $tumor \
                --out_file $out_file


done < $TNFILE


## Controls
while read tumor normal; do
    
    tn=$tumor--$normal
    plasma_vcf=$WORKING_DIR/control/$tn.union.vcf
    out_file=$WORKING_DIR/summary/control/$tumor.detection_summary.txt

    Rscript $SRCDIR/init-detection-summary.r \
        --clones $CLONE_ASSIGNMENTS \
        --vcf $plasma_vcf \
        --sample_id $tumor \
        --out_file $out_file

done < $CONTROLS



#########################
## Downstream analysis ##
#########################

## Merge, compute noise floor, estimate tumor fraction
Rscript $SRCDIR/init-tf-estimation.r \
    --in_dir $WORKING_DIR/summary \
    --control_dir $WORKING_DIR/summary/control \
    --id_map $ID_MAP \
    --tn_file $TNFILE \
    --out_file $WORKING_DIR/WCM63-tumor-fraction-estimates.txt


## Integrate tumor copy number and mutation multiplcity information 
Rscript $SRCDIR/init-mean-cn-multiplicity.r \
    --clonal_clone_ids 0,1 \
    --in_file_pyclone $PYCLONE_DIR/WCM63.pyclone_vi_input.conipher_converted.outliers_excluded.txt \
    --in_file_phyclone $CLONE_ASSIGNMENTS \
    --in_dir_pileup $WORKING_DIR/pileup \
    --in_file_detection_rate $WORKING_DIR/WCM63-tumor-fraction-estimates.txt \
    --pn_file $TNFILE \
    --out_file $WORKING_DIR/WCM63-eclipse-input.txt


## Compute a cn/multiplicity-adjusted tumor fraction for each plasma 
## sample using truncal clones
Rscript $SRCDIR/init-clonal-adjusted-tf-estimation.r \
    --in_file_det_rate $WORKING_DIR/WCM63-tumor-fraction-estimates.txt \
    --in_file_cn $WORKING_DIR/WCM63-eclipse-input.txt \
    --out_file $WORKING_DIR/WCM63-adjusted-clonal-tf-estimates.txt


## Compare cosine similarity between tumor/plasma read depth profiles
Rscript $SRCDIR/init-compare-ichor-fragcounter.r \
    --in_dir_cfdna $CFDNA_DIR \
    --in_dir_tumor $PROJECT_DIR/jabba-latest/dryclean/tumor-decomposition \
    --tn_file $TNFILE_TUMOR \
    --pn_file $TNFILE \
    --out_file $WORKING_DIR/ichor-fragcounter-read-depth-comparison.txt 
