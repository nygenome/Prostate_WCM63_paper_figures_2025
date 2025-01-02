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
## Take the union of jabba event footprints of a specific type, overlapping a specific chromosome
libs = c('optparse', 'GenomicRanges', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Convert a footprint record into GRanges
fp2gr = function(x, id) {
  
  x = gsub('\\+(?=,)|\\-(?=,)','',x, perl=T)
  x = gsub('\\+(?=;)|\\-(?=;)','',x, perl=T)
  x = gsub('\\+$|\\-$','',x, perl=T)
  x = strsplit(x,';|,')
  
  res = list()
  
  for (i in 1:length(x)) {
    
    seqnames = gsub(':.*','',x[[i]])
    start = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[1]))
    end = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[2]))
    
    res[[i]] = GRanges(seqnames, IRanges(start=start, end=end))
  }
  
  res = do.call(c, res)
  res$id = id
  
  return(res) 
  
}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir"),     type='character', help="Directory holding JaBbA results"),
  make_option(c("-t", "--tn_file"),    type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-c", "--chr"),        type='character', help="Comma-delimited list of bed files to select footprints"),
  make_option(c("-f", "--force_chr"),  type='logical',   help="When TRUE, only include intervals overlapping --chr", default=FALSE),
  make_option(c("-e", "--event_type"), type='character', help="Comma-delimited list of tumor IDs (will automatically pull normals as controls)"),
  make_option(c("-o", "--out_file"),   type='character', help="Output file"))
opt = parse_args(OptionParser(option_list=option_list))
opt$event_type = unlist(strsplit(opt$event_type, ',', fixed=T))


## Read tumor normal pairs, get input files
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor', 'normal'))
tn$pair_name = paste0(tn$tumor, '--', tn$normal)
tn$in_file = paste0(opt$in_dir, '/', tn$pair_name, '/jabba.events.genes.tab')

if(!all(file.exists(tn$in_file))) stop('Missing files!')


## Read footprint tables
fp = lapply(tn$in_file, read.csv, h=T, stringsAsFactors=F, sep='\t')


## Select event of interest
fp = lapply(fp, function(x) x[x$type %in% opt$event_type & grepl(paste0(opt$chr,':'),x$footprint, fixed=T), c('type', 'footprint')])
fp = do.call(rbind, fp)


## Convert to GRanges
res = Reduce(union, mapply(FUN=fp2gr, x=fp$footprint, id=fp$type))
res = as.data.frame(res)[, c('seqnames', 'start', 'end')]


## Filter if specified 
if (opt$force_chr) {
  res = res[res$seqnames %in% opt$chr, ]
}


write.table(res, opt$out_file, row.names=F, col.names=F, quote=F, sep='\t')
message(opt$out_file)
