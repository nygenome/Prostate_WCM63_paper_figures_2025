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
## Integrate multiplicity/cn data
libs = c('optparse', 'data.table', 'reshape2', 'VariantAnnotation')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



read.vcf = function(f, sample_id) {

  x = readVcf(f)
  res = rowRanges(x)
  res$alt_reads = sapply(geno(x)$AD[, sample_id], `[`, 2)
  res$total_reads = geno(x)$DP[, sample_id]

  seqlevelsStyle(res) = 'NCBI'

  names(res) = paste0(as.character(seqnames(res)), ':', start(res))

  return(res)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file_pyclone"),        type='character', help="Pyclone input file"),
  make_option(c("-I", "--in_file_phyclone"),       type='character', help="TSV output from Phyclone"),
  make_option(c("-p", "--in_dir_pileup"),          type='character', help="Plasma pileup edir"),
  make_option(c("-d", "--in_file_detection_rate"), type='character', help="Sample ID"),
  make_option(c("-t", "--pn_file"),                type='character', help="Plasma-normal pairs"),
  make_option(c("-c", "--clonal_clone_ids"),       type='character', help="Comma-delimited list of clonal clones (e.g. present in all samples)"),
  make_option(c("-o", "--out_file"),               type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))



#####################
## Read input data ##
#####################

opt$clonal_clone_ids = unlist(strsplit(opt$clonal_clone_ids, ',', fixed=T))

pn = fread(opt$pn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('plasma', 'normal'), data.table=F)
pn$pair_name = paste0(pn$plasma, '--', pn$normal)
pn$in_file = paste0(opt$in_dir_pileup,'/',pn$pair_name,'.union.vcf')

det = fread(opt$in_file_detection_rate, h=T, stringsAsFactors=F, sep='\t', data.table=F)
det = det[!duplicated(det$clone), c('clone', 'mean_error_rate')]
pyc = fread(opt$in_file_pyclone, h=T, stringsAsFactors=F, sep='\t', data.table=F)
phy = fread(opt$in_file_phyclone, h=T, stringsAsFactors=F, sep='\t', data.table=F)


## Merge pyclone outputs
phy = merge(x=pyc, y=phy, on=c('mutation_id', 'sample_id'))

## Clone to mutation mapping 
cln = phy[!duplicated(phy$mutation_id), c('mutation_id', 'clone_id')]



##########################
## Compute multiplicity ##
##########################

phy$total_reads = phy$alt_counts + phy$ref_counts
phy$vaf = phy$alt_counts / phy$total_reads

phy$total_cn = phy$major_cn + phy$minor_cn
phy$avg_cn = (phy$tumour_content * phy$total_cn) + ((1 - phy$tumour_content) * phy$normal_cn)

## Disallow multiplicities <0 or >maj cn
phy$multiplicity = round((phy$vaf * phy$avg_cn) / phy$tumour_content)
phy$multiplicity = pmax(phy$multiplicity, 1)
phy$multiplicity = pmin(phy$multiplicity, phy$major_cn)

## Average purity, vaf, ccf, copy number and multiplicity across all tumor samples that 
## have a given mutation
phy$tumour_mean_ccf = phy$ccf
phy$tumour_mean_vaf = phy$vaf
phy = phy[, c('mutation_id', 'sample_id', 'normal_cn', 'total_cn', 'multiplicity', 'tumour_content', 'tumour_mean_ccf', 'tumour_mean_vaf')]
phy = melt(phy, id.var=c('mutation_id', 'sample_id'))
phy = dcast(mutation_id ~ variable, data=phy, value.var='value', fun.aggregate=mean, na.rm=T)



#########################
## Clone-specific info ##
#########################

## Add clone-specific error rate
phy$clone = cln$clone_id[match(phy$mutation_id, cln$mutation_id)]
phy$error_rate = det$mean_error_rate[match(phy$clone, det$clone)]

## Mark "clonal" clones, e.g. those present across all samples
phy$is_clonal = phy$clone %in% opt$clonal_clone_ids

phy$mutation_id_short = gsub('\\:[[:upper:]].*', '', phy$mutation_id)
phy[, c('chr', 'pos', 'ref', 'alt')] = colsplit(phy$mutation_id, pattern=':', names=c('chr', 'pos', 'ref', 'alt'))



######################
## Read plasma VCFs ##
######################

res = c()
for (i in 1:nrow(pn)) {

  message(pn$plasma[i])

  vcf = read.vcf(f=pn$in_file[i], sample_id=pn$plasma[i])

  ## Annotate plasma with tumor/cluster-informed info from above 
  ## 
  ## Some clusters were excluded from tree building, 
  ## and get_vaf.py is unable to compute read evidence 
  ## for long indels so this number won't line up exactly
  res.i = phy
  res.i$sample_id = pn$plasma[i]
  res.i$alt_reads = vcf$alt_reads[match(res.i$mutation_id_short, names(vcf))]
  res.i$total_reads = vcf$total_reads[match(res.i$mutation_id_short, names(vcf))]
  res.i = res.i[!is.na(res.i$alt_reads) & !is.na(res.i$alt_reads), ]

  res = rbind(res, res.i)

}


col.order = c('sample_id', 'mutation_id', 'chr', 'pos', 'ref', 'alt', 'alt_reads', 
              'total_reads', 'total_cn', 'multiplicity', 'normal_cn', 'clone', 
              'tumour_content', 'tumour_mean_ccf', 'tumour_mean_vaf', 'error_rate', 
              'is_clonal')
res = res[, col.order]


## Write result 
write.table(res, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
