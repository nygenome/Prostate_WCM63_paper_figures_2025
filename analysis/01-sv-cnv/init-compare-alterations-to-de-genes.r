#!/nfs/sw/R/R-4.0.0/bin/Rscript
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
## Compare WGS alterations to DE genes 
libs = c('optparse', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Tally alterations by histology. histology should be aligned with 
## columns of x
tally.by.histology = function(x, histology) {

  uniq.hist = unique(histology)
  res = as.data.frame(matrix(NA, ncol=length(uniq.hist)+1, nrow=nrow(x)))
  colnames(res) = c('gene', uniq.hist)
  res$gene = rownames(x)

  ## Tally the alterations in each histology
  for (i in 1:length(uniq.hist)) {
      res[, uniq.hist[i]] = apply(x[, histology == uniq.hist[i]], 1, function(y) sum(!is.na(y)))
  }

  return(res)

}



## Summarize alterations in relation to a table of DE genes
summarize.by.histology = function(x, de, histology.counts) {

  ## Conserved in a histology and absent from the other 
  for (i in names(histology.counts)) {

    absent.in.others =  x[, setdiff(names(histology.counts), i)] == 0

    ## Strict definition of conservation
    conserved.i = x[, i] == histology.counts[i]
    message('Fully conserved in ', i, ' and absent in others: ', sum(conserved.i & absent.in.others))

    if (sum(conserved.i & absent.in.others) > 0) {
      message('Overlap with DE gene list:')
      overlap = intersect(rownames(de), x$gene[conserved.i & absent.in.others])
      if (length(overlap) > 0) {
        res = cbind(x[match(overlap, x$gene), ], de[overlap, c('log2FoldChange', 'padj')])
        write.table(res, row.names=F, col.names=T, quote=F, sep='\t')
      }

    }

    ## Relaxed definition
    min.samples = floor((2/3) * histology.counts[i])
    conserved.i = x[, i] >= min.samples
    message('Partially conserved in ', i, ' (>=', min.samples, ' samples) and absent in others: ', sum(conserved.i & absent.in.others))

    if (sum(conserved.i & absent.in.others) > 0) {
      message('Overlap with DE gene list:')

      overlap = intersect(rownames(de), x$gene[conserved.i & absent.in.others])
       if (length(overlap) > 0) {
        res = cbind(x[match(overlap, x$gene), ], de[overlap, c('log2FoldChange', 'padj')])
        write.table(res, row.names=F, col.names=T, quote=F, sep='\t')
      }
    }

    message('\n')
    
  } 

}



## Get arguments
option_list = list(
  make_option(c("-i", "--onco_mtx"), type='character', help="Oncoprint app input matrix"),
  make_option(c("-t", "--de_genes"), type='character', help="DE genes"),
  make_option(c("-m", "--metadata"), type='character', help="Sample metadata"))
opt = parse_args(OptionParser(option_list=option_list))



FDR_CUTOFF = 0.01


## Read data 
onco = read.csv(opt$onco_mtx, h=T, stringsAsFactors=F, sep='\t')
onco$Impact[onco$Impact == ''] = NA
onco$CNV_Class[onco$CNV_Class == ''] = NA
onco$CNV_Class[onco$CNV_Scale == 'largescale'] = NA
onco$CNV_Class[onco$CNV_Class == 'deletion'] = 'loss'
onco$CNV_Class[onco$CNV_Class == 'amplification'] = 'gain'

mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta$pair_name = paste0(mta$tumor, '--', mta$normal)

de = read.csv(opt$de_genes, h=T, stringsAsFactors=F, row.names=1)
de = de[!is.na(de$padj) & de$padj < FDR_CUTOFF, ]



#######################
## Compare SNV/INDEL ##
#######################

snv = dcast(formula=Gene ~ Sample, value.var='Impact', data=onco)
rownames(snv) = snv$Gene
snv = snv[, !colnames(snv) %in% 'Gene']

message('\n\n-----SNV/INDEL-----')
histology = mta$clinical_classification[match(colnames(snv), mta$pair_name)]
snv = tally.by.histology(x=snv, histology=histology)
summarize.by.histology(x=snv, de=de, histology.counts=table(histology))
message('\n\n')



###################
## Compare gains ##
###################

gain = onco
gain$CNV_Class[gain$CNV_Class != 'gain'] = NA
gain = dcast(formula=Gene ~ Sample, value.var='CNV_Class', data=gain)
rownames(gain) = gain$Gene
gain = gain[, !colnames(gain) %in% 'Gene']

message('\n\n-----CNV: gain-----')
histology = mta$clinical_classification[match(colnames(gain), mta$pair_name)]
gain = tally.by.histology(x=gain, histology=histology)
summarize.by.histology(x=gain, de=de, histology.counts=table(histology))
message('\n\n')



####################
## Compare losses ##
####################

loss = onco
loss$CNV_Class[loss$CNV_Class != 'loss'] = NA
loss = dcast(formula=Gene ~ Sample, value.var='CNV_Class', data=loss)
rownames(loss) = loss$Gene
loss = loss[, !colnames(loss) %in% 'Gene']

message('\n\n-----CNV: loss-----')
histology = mta$clinical_classification[match(colnames(loss), mta$pair_name)]
loss = tally.by.histology(x=loss, histology=histology)
summarize.by.histology(x=loss, de=de, histology.counts=table(histology))
message('\n\n')
