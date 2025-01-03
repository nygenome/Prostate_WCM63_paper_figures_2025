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
## Merge detection rates for plasma and unrelated (or healthy) controls
## Compare detection rate against controls to determine whether plasma signal 
## is above lower limit of detection
libs = c('optparse', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir"),      type='character', help="TSV output from Phyclone"),
  make_option(c("-c", "--control_dir"), type='character', help="VCF output from get_vaf.py"),
  make_option(c("-t", "--tn_file"),     type='character', help="Sample ID"),
  make_option(c("-m", "--id_map"),      type='character', help="ID map"),
  make_option(c("-o", "--out_file"),    type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))


Z_SCORE_CUTOFF = 1.2


## Read tumor normal pairs 
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal'))
tn$in_file = paste0(opt$in_dir, '/', tn$tumor, '.detection_summary.txt')


## ID mapping 
id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')


## Read input tables 
dta = do.call(rbind, lapply(tn$in_file, read.csv, h=T, stringsAsFactors=F, sep='\t'))
head(dta)
dta$internal_sample_id = dta$sample_id
dta$sample_id = id.map$New.sample.name[match(dta$sample_id, id.map$Associated.sample)]


## Read control data
ctl.files = Sys.glob(paste0(opt$control_dir,'/*.detection_summary.txt'))
ctl = do.call(rbind, lapply(ctl.files, read.csv, h=T, stringsAsFactors=F, sep='\t'))


## Compute sample-wide and clone-specific error rates 
error.rate = as.data.frame(do.call(rbind, tapply(ctl$detection_rate, ctl$clone, function(x) c(mean(x), sd(x)))))
colnames(error.rate) = c('mean_error_rate', 'sd_error_rate')


## Handle error rate of 0 (no variant supporting reads detected)
if (any(error.rate$mean_error_rate == 0)) {
  error.rate[which(error.rate$mean_error_rate == 0), ] = error.rate['all', ]
}


## Compute z-scores 
dta$z_score = (dta$detection_rate - dta$mean_error_rate) / dta$sd_error_rate
dta$z_score_cutoff = Z_SCORE_CUTOFF
dta$detected = dta$z_score > Z_SCORE_CUTOFF

ctl$z_score = (ctl$detection_rate - ctl$mean_error_rate) / ctl$sd_error_rate
ctl$detected = ctl$z_score > Z_SCORE_CUTOFF


## Map z-score cutoff to TF LLOD 
dta$tf_cutoff = (dta$z_score_cutoff * dta$sd_error_rate) + dta$mean_error_rate


## Write result 
write.table(dta, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
