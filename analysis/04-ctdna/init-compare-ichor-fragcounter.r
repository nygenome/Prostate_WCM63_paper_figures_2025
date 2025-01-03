#!/nfs/sw/R/R-4.0.0/bin/Rscript
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# William F. Hooper

################################################################# /COPYRIGHT ###
################################################################################
## Compare cosine similarity between read depth profiles
libs = c('optparse', 'data.table', 'GenomicRanges', 'gUtils', 'lsa')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


read.ichor = function(f) { 

  x = fread(f, h=T, stringsAsFactors=F, sep='\t', data.table=F)
  x = makeGRangesFromDataFrame(x, keep.extra=T)

  colnames(mcols(x)) = gsub('^.*\\.', '', colnames(mcols(x)))

  return(x)

}



option_list = list(
  make_option(c("-c", "--in_dir_cfdna"),    type='character', help="Input directory"),
  make_option(c("-t", "--in_dir_tumor"),    type='character', help="Input directory"),
  make_option(c("-v", "--tn_file"),         type='character', help="Tumor normal pairs"),
  make_option(c("-s", "--pn_file"),         type='character', help="Plasma normal pairs"),
  make_option(c("-o", "--out_file"),        type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))



## Read pairs 
tn = fread(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor', 'normal'), data.table=F)
pn = fread(opt$pn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('plasma', 'normal'), data.table=F)

tn$in_file = paste0(opt$in_dir_tumor,'/',tn$tumor,'.dryclean-decomp.rds')
pn$in_file = paste0(opt$in_dir_cfdna,'/Sample_',pn$plasma,'/compbio/ichorCNA_default/',pn$plasma,'.cna.seg')



## Read data 
tumor.cnv = lapply(tn$in_file, readRDS)
names(tumor.cnv) = tn$tumor 

plasma.cnv = lapply(pn$in_file, read.ichor)
names(plasma.cnv) = pn$plasma 


## Compute cosine similarity between all tumor/plasma combinations
res = matrix(NA, 
             nrow=length(tumor.cnv), 
             ncol=length(plasma.cnv), 
             dimnames=list(names(tumor.cnv), names(plasma.cnv)))

for (i in 1:nrow(res)) {
  for (j in 1:ncol(res)) {

    cnv = plasma.cnv[[j]] %$% tumor.cnv[[i]]
    cnv = cnv[!is.na(cnv$foreground.log) & !is.na(cnv$logR)]

    res[i, j] = cosine(x=cnv$foreground.log, y=cnv$logR)

  }
}


## Write result
write.table(res, opt$out_file, row.names=T, col.names=T, quote=F, sep='\t')
message(opt$out_file)
