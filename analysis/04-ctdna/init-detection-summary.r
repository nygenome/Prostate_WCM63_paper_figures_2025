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
## Compute detection rates for a set of tumor-derived clones 
libs = c('optparse', 'data.table', 'reshape2', 'VariantAnnotation', 'boot')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



read.vcf = function(f, sample_id) {

  x = readVcf(f)
  res = rowRanges(x)
  res$alt_reads = sapply(geno(x)$AD[, sample_id], `[`, 2)
  res$total_reads = geno(x)$DP[, sample_id]

  names(res) = paste0(gsub('^chr', '', as.character(seqnames(res))), ':', start(res))

  return(res)

}



detection.rate = function(x, i) {

  return(sum(x$alt_reads[i]) / sum(x$total_reads[i]))

}



## Get arguments
option_list = list(
  make_option(c("-c", "--clones"),    type='character', help="TSV output from Phyclone"),
  make_option(c("-v", "--vcf"),       type='character', help="VCF output from get_vaf.py"),
  make_option(c("-s", "--sample_id"), type='character', help="Sample ID"),
  make_option(c("-o", "--out_file"),  type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))



## Read mutation --> clone assignments
cln = fread(opt$clones, h=T, stringsAsFactors=F, sep='\t', data.table=F)
cln = cln[!duplicated(cln$mutation_id), c('mutation_id', 'clone_id'), ]
cln[, c('chr', 'pos', 'ref', 'alt')] = colsplit(cln$mutation_id, pattern=':', names = c('chr', 'pos', 'ref', 'alt'))
cln$mutation_id = paste0(as.character(cln$chr),':',cln$pos)


## Read plasma mutations
vcf = read.vcf(opt$vcf, opt$sample_id)


## Harmonize mutation list
## Some clusters were excluded from tree building, 
## and get_vaf.py is unable to compute read evidence 
## for long indels 
mut.keep = intersect(cln$mutation_id, names(vcf))
cln = cln[cln$mutation_id %in% mut.keep, ]
vcf = vcf[mut.keep]
message('Kept ', length(mut.keep), ' mutations')


## Annotate mutations with cluster membership 
vcf$clone_id = cln$clone_id[match(names(vcf), cln$mutation_id)]


res = data.frame(sample_id=opt$sample_id,
                 clone=c('all', unique(vcf$clone_id)),
                 num_variants=NA,
                 alt_reads=NA,
                 reads_checked=NA,
                 mean_depth=NA,
                 detection_rate=NA)

for (i in 1:nrow(res)) {

  clone = res$clone[i]

  if (clone == 'all') {
    vcf.clone = vcf
  } else {
    vcf.clone = vcf[vcf$clone_id == clone]
  }

  ## Compute a 95% confidence interval for the detection rate 
  clone.boot = boot(data=vcf.clone, 
                    statistic=detection.rate, 
                    R=1000)

  ci = boot.ci(boot.out=clone.boot,
               conf=0.95, 
               type='norm')


  res$num_variants[i] = length(vcf.clone)
  res$alt_reads[i] = sum(vcf.clone$alt_reads)
  res$reads_checked[i] = sum(vcf.clone$total_reads)
  res$mean_depth[i] = mean(vcf.clone$total_reads)
  res$detection_rate[i] = res$alt_reads[i] / res$reads_checked[i]
  res$detection_rate.ci95.lower[i] = ci$normal[, 2]
  res$detection_rate.ci95.upper[i] = ci$normal[, 3]

}


## Lower CI shouldn't be below 0 
res$detection_rate.ci95.lower[res$detection_rate.ci95.lower < 0] = 0


## Write result
write.table(res, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
