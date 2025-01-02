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
## Reformat mutation/copy number data for input into PyClone-VI
libs = c('optparse', 'GenomicRanges', 'gUtils', 'VariantAnnotation')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200)



## Get arguments
option_list = list(
  make_option(c("-t", "--tn_file"),       type='character', help="Tumor normal pairs file"),
  make_option(c("-c", "--cn_dir"),        type='character', help="Input jabba.simple.rds with slot agtrack (i.e., jabba was run with a het pileup supplied"),
  make_option(c("-p", "--purity_ploidy"), type='character', help="Tumor purity/ploidy table"),
  make_option(c("-v", "--vcf_dir"),       type='character', help="Directory with VCFs"),
  make_option(c("-g", "--gender"),        type='character', help="Gender"),
  make_option(c("-o", "--out_file"),      type='character', help="Output RDS with per-segment allele-specific copy number"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor normal pairs 
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

## Read purity/ploidy
pp = read.csv(opt$purity_ploidy, h=F, stringsAsFactors=F, col.names=c('tumor', 'purity', 'ploidy'))
tn$purity = pp$purity[match(tn$tumor, pp$tumor)]

## Get input files 
tn$cn_file = paste0(opt$cn_dir,'/',tn$pair_name,'-jabba-ascn.rds')
tn$vcf = paste0(opt$vcf_dir,'/',tn$pair_name,'.union.vcf')



if(!all(file.exists(tn$cn_file))) {
  stop("Missing copy number data")
}
if(!all(file.exists(tn$vcf))) {
  stop("Missing somatic VCFs")
}



col.sel = c('mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 
            'normal_cn', 'major_cn', 'minor_cn', 'tumour_content')
res = c()

for (i in 1:nrow(tn)) {
  message(tn$tumor[i])

  ## Read data, filter VCF for high conf / union variants
  cn = readRDS(tn$cn_file[i])
  vcf = readVcf(tn$vcf[i])
  vcf = vcf[info(vcf)$HighConfidence | info(vcf)$TYPE == 'ADDED']

  ## Extract ref/alt counts
  rowRanges(vcf)$ref_counts = sapply(geno(vcf)$AD[, tn$tumor[i]], `[`, 1)
  rowRanges(vcf)$alt_counts = sapply(geno(vcf)$AD[, tn$tumor[i]], `[`, 2)

  ## Add normal cn 
  rowRanges(vcf)$normal_cn = 2
  if (opt$gender == 'male') {
    rowRanges(vcf)$normal_cn[as.character(seqnames(rowRanges(vcf))) %in% c('chrX', 'chrY', 'X', 'Y')] = 1
  } 

  ## Add major/minor cn 
  vcf = rowRanges(vcf)
  vcf = vcf[!duplicated(names(vcf))]
  vcf = as.data.frame(gr.val(vcf, cn, c('major_cn', 'minor_cn')))

  ## Add some additional columns 
  vcf$mutation_id = rownames(vcf)
  vcf$tumour_content = tn$purity[i]
  vcf$sample_id = tn$tumor[i]
  rownames(vcf) = c()

  vcf = vcf[, col.sel]
  res = rbind(res, vcf)

}


## Exclude any mutations that are missing copy number data 
## in any samples 
mut.exclude = res$mutation_id[is.na(res$major_cn) | is.na(res$minor_cn)]
if (length(mut.exclude) > 0) {
  res = res[!res$mutation_id %in% mut.exclude, ]
}


write.table(res, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
