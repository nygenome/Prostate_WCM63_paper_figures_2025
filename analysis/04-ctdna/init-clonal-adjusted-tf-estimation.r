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
## Compute a multiplicity/cn adjusted tumor fraction using eqn b in  
## https://www.nature.com/articles/s41586-023-05776-4/figures/12
libs = c('optparse', 'data.table')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)




## Get arguments
option_list = list(
  make_option(c("-i", "--in_file_det_rate"), type='character', help="TF estimates"),
  make_option(c("-t", "--in_file_cn"),       type='character', help="Data formatted for ECLIPSE"),
  make_option(c("-o", "--out_file"),         type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data
det = fread(opt$in_file_det_rate, h=T, stringsAsFactors=F, sep='\t', data.table=F)
cn = fread(opt$in_file_cn, h=T, stringsAsFactors=F, sep='\t', data.table=F)


## Subset to clone of interest
det = det[det$clone == '0', ]
cn = cn[cn$clone == '0', ]


## Get mean normal cn, total cn, multiplicity 
det$mean_ncn = mean(cn$normal_cn)
det$mean_tcn = mean(cn$total_cn)
det$mean_mcn = mean(cn$multiplicity)
det$adjusted_tf = det$mean_ncn / ((det$mean_mcn / det$detection_rate) - det$mean_tcn + det$mean_ncn)


## Write result
res = det[, c('sample_id', 'clone', 'detection_rate', 'adjusted_tf')]
write.table(res, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)