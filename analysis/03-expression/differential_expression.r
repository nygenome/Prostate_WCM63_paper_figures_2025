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

library(DESeq2)

#Read in batch-corrected but NOT library-size-normalized counts.
#Since we used ComBat_seq, these results can be used as input to DESeq2.

#Looks like counts including WCM0 were used as input to DE, but counts excluding WCM0 were used for expression heatmap.
#Difference should be very minor, though.

exp <- read.csv("Results_batch_correction/NEPC_GRCh37_counts_corrected_incl_PM0_as_polyA.csv",header=TRUE,row.names=1,check.names=FALSE)
exp <- exp[,grep('WCM63',colnames(exp))]

#Set up design table specifying whether each sample is squamous or adeno.

WCM63_sample_ids <- c("WCM63_B","WCM63_A","WCM63_E","WCM63_F","WCM63_G","WCM63_I","WCM63_J","WCM63_K","WCM63_L","WCM63_D","WCM63_H","WCM63_C")
WCM63_groups <- c("SCC","Ad","Ad","Ad","Ad","SCC","SCC","SCC","SCC","Ad","SCC","Ad")

design <- data.frame(sample = WCM63_sample_ids,
    group = WCM63_groups,
    row.names = WCM63_sample_ids,
    check.names=FALSE)

design <- design[colnames(exp),]

design <- data.frame(group = design$group,
    row.names=design$sample,
    check.names=FALSE)

#Run differential expression using a simple pairwise design.

dds <- DESeq(DESeqDataSetFromMatrix(countData=exp,colData=design,design=~group))
results_squam_vs_adeno <- results(dds,contrast=c("group","Squamous","Adeno"))
results_squam_vs_adeno <- as.data.frame(results_squam_vs_adeno)

#Print to CSV.

write.csv(results_squam_vs_adeno,file="DE_PM63_squam_vs_adeno_after_batch_correct.csv")

#Let's also print the "rnk" file needed as input to GSEA.
#We will use the "stat" column, which basically incorporates information about the magnitude of change, the direction of change (+/- log-fold-change), and the significance.
#We also filter out very low expression values (require baseMean >= 30) and NA p-values.

rnk_ind <- which(!is.na(results_squam_vs_adeno$padj) & results_squam_vs_adeno$baseMean >= 30)
rnk <- results_squam_vs_adeno[rnk_ind,]
rnk <- rnk[order(rnk$stat),]
rnk <- data.frame(gene = rownames(rnk),stat = rnk$stat)

write.table(rnk,file="DE_PM63_squam_vs_adeno_after_batch_correct.rnk",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
