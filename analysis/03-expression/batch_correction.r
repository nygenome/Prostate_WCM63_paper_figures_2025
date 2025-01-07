#!/nfs/sw/R/R-4.0.0/bin/Rscript
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

#Load libraries.

library(sva)
library(DESeq2)

#Read in counts for the full rapid autopsy cohort, which includes WCM63 along with several other patients.
#Also included WCM12 as well - not sure why this is noted as a separate thing from just the RA cohort. Getting confirmation on this now.

counts <- read.csv("raw_counts_RA_and_PM12_only.csv",header=TRUE,row.names=1,check.names=FALSE)

#Set up design table for WCM63 with prep information (polyA-selection vs. ribo-depletion).

WCM63_sample_ids <- c("WCM63_B","WCM63_A","WCM63_E","WCM63_F","WCM63_G","WCM63_I","WCM63_J","WCM63_K","WCM63_L","WCM63_D","WCM63_H","WCM63_C")
WCM63_prep_for_sample_ids <- c(rep("polyA",times=6),
    rep("ribodepletion",times=5),
    "polyA")

design_WCM63 <- data.frame(sample = WCM63_sample_ids,
    prep = WCM63_prep_for_sample_ids,
    row.names=WCM63_sample_ids,
    check.names=FALSE)

#Read in a csv file with prep information for each of the non-WCM63 samples.
#Here is an example of what this file looks like for a subset of a few patients (sample IDs are placeholders and do not necessarily correspond to the real IDs):

# "sample","prep"
# "WCM12_A","polyA"
# "WCM12_B","polyA"
# "WCM159_A","ribodepletion"
# "WCM159_B","ribodepletion"
# "WCM159_C","ribodepletion"
# "WCM159_D","ribodepletion"
# "WCM677_A","polyA"
# "WCM677_B","ribodepletion"
# "WCM677_C","ribodepletion"
# "WCM677_D","ribodepletion"
# "WCM677_E","ribodepletion"
# "WCM677_F","ribodepletion"
# "WCM1070_A","polyA"
# "WCM1070_B","polyA"
# "WCM1070_C","polyA"
# "WCM1070_D","polyA"
# "WCM1070_E","polyA"
# "WCM1089_A","polyA"
# "WCM1089_B","polyA"
# "WCM1089_C","polyA"
# "WCM1089_D","ribodepletion"

design_non_WCM63 <- read.csv("prep_non_PM63.csv",header=TRUE,check.names=FALSE)
rownames(design_non_WCM63) <- design_non_WCM63$sample

#The only patient that this file does not contain information on is WCM0.
#We were never able to 100% confirm the prep information for samples from this patient with the lab.
#However, the expression pattern and clustering of samples from this patient strongly suggests that they were all done with polyA prep.
#So, we try running batch correction both including and excluding this patient, to see how the results compare.

design_non_WCM63_excl_WCM0 <- design_non_WCM63

WCM0_prep_to_rbind <- data.frame(sample = c("WCM0_A","WCM0_B","WCM0_C","WCM0_D","WCM0_E","WCM0_F"),
    prep = "polyA",
    row.names=c("WCM0_A","WCM0_B","WCM0_C","WCM0_D","WCM0_E","WCM0_F"),
    check.names=FALSE)

design_non_WCM63_incl_WCM0 <- rbind(design_non_WCM63,WCM0_prep_to_rbind)

#Let's also make a version of the counts table both with and without WCM0.
#And combine prep info for WCM63 and non-WCM63 samples, and match up for sample IDs to be in the same order in both the counts and the design table.

counts_incl_WCM0 <- counts
counts_excl_WCM0 <- counts[,grep('WCM0',colnames(counts),invert=TRUE)]

design_incl_WCM0 <- rbind(design_WCM63,design_non_WCM63_incl_WCM0)
design_excl_WCM0 <- rbind(design_WCM63,design_non_WCM63_excl_WCM0)

design_incl_WCM0 <- design_incl_WCM0[colnames(counts_incl_WCM0),]
design_excl_WCM0 <- design_excl_WCM0[colnames(counts_excl_WCM0),]

#Ready to run batch correction.

counts_corrected_minus_WCM0 <- ComBat_seq(counts=as.matrix(counts_excl_WCM0),batch=design_excl_WCM0$prep)
counts_corrected_incl_WCM0_as_polyA <-  ComBat_seq(counts=as.matrix(counts_incl_WCM0),batch=design_incl_WCM0$prep)

write.csv(counts_corrected_minus_WCM0,file="Results_batch_correction/NEPC_GRCh37_counts_corrected_minus_PM0.csv",quote=FALSE)
write.csv(counts_corrected_incl_WCM0_as_polyA,file="Results_batch_correction/NEPC_GRCh37_counts_corrected_incl_PM0_as_polyA.csv",quote=FALSE)

#Subset to WCM63 only from each of these matrices, then run library size normalization using DESeq2.

counts <- counts_corrected_minus_WCM0[,grep('WCM63',colnames(counts_corrected_minus_WCM0))]
samples <- colnames(counts)
design <- data.frame(s = samples,row.names=samples,check.names=FALSE)
cds <- DESeqDataSetFromMatrix(countData=counts,colData=design,design=~1)
cds <- estimateSizeFactors(cds)
norm_counts <- counts(cds,normalized=TRUE)

write.csv(norm_counts,file="Results_batch_correction/NEPC_PM63_GRCh37_counts_libnorm_batch_corrected_using_RA_cohort_minus_PM0.csv")

counts <- counts_corrected_incl_WCM0_as_polyA[,grep('WCM63',colnames(counts_corrected_incl_WCM0_as_polyA))]
samples <- colnames(counts)
design <- data.frame(s = samples,row.names=samples,check.names=FALSE)
cds <- DESeqDataSetFromMatrix(countData=counts,colData=design,design=~1)
cds <- estimateSizeFactors(cds)
norm_counts <- counts(cds,normalized=TRUE)

write.csv(norm_counts,file="Results_batch_correction/NEPC_PM63_GRCh37_counts_libnorm_batch_corrected_using_RA_cohort_incl_PM0.csv")
