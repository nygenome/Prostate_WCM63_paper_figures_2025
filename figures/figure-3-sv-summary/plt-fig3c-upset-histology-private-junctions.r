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
## Generate upset plots comparing histology-private junctions
libs = c('optparse', 'UpSetR', 'grid')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



option_list = list(
  make_option(c("-i", "--in_file"),   type='character', help="Output of analysis/01-sv-cnv/patient-junction-support-by-bed.r"),
  make_option(c("-t", "--tn_file"),   type='character', help="Tumor-normal pairs file"),
  make_option(c("-m", "--metadata"),  type='character', help="Sample metadata"),
  make_option(c("-M", "--id_map"),    type='character', help="Map to manuscript IDs (optional)"),
  make_option(c("-o", "--out_file"),  type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor-normal pairs
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t')

if (ncol(tn) == 2) {
  colnames(tn) = c('tumor', 'normal')
} else {
  colnames(tn) = c('tumor','normal','gender')
}
tn$pair_name = paste0(tn$tumor,'--',tn$normal)


## Optionally update names
if (!is.null(opt$id_map)) {

    id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
    tn$name = id.map$New.sample.name[match(gsub('-','_',tn$tumor), id.map$Associated.sample)]

} else {

  tn$name = tn$tumor

}


## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta$pair_name = paste0(mta$tumor,'--',mta$normal)
tn[, c('site', 'histology')] = mta[match(tn$pair_name, mta$pair_name), c('simplified_site', 'clinical_classification')]


## Read data and align with TN pairs
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, check.names=F, row.names=1, sep='\t')
dta = dta[, tn$tumor]
names(dta) = tn$name


## Convert to binary matrix, anything with at least one 
## read is considered "present"
for (i in 1:ncol(dta)) {

  dta[, i] = as.numeric(dta[, i] > 0)

}


## Identify squam- and adeno-private junctions
squam.private = apply(dta[, tn$histology == 'CRPC-Ad'], 1, function(x) all(x == 0)) ## Private to squam, i.e. no read support in all adenos
adeno.private = apply(dta[, tn$histology == 'CRPC-SCC'], 1, function(x) all(x == 0)) ## Private to adeno, i.e. no read support in all squams

pct.squam.private = paste0(round(sum(squam.private) / nrow(dta), 4) * 100,'%')
pct.squam.private = paste0('CRPC-SCC-private junctions\n', pct.squam.private, ' (',sum(squam.private),'/',nrow(dta),')')

pct.adeno.private = paste0(round(sum(adeno.private) / nrow(dta), 4) * 100,'%')
pct.adeno.private = paste0('CRPC-Ad-private junctions\n',pct.adeno.private, ' (',sum(adeno.private),'/',nrow(dta),')')



##########
## Plot ##
##########

pdf(opt$out_file)

## Squamous-private
upset(data=dta[squam.private, tn$histology == 'CRPC-SCC'], 
      order.by='freq', 
      nsets=sum(tn$histology == 'CRPC-SCC'),
      main.bar.color='#CC61B0',
      sets.bar.color='#CC61B0',
      text.scale=1.9)
grid.text(pct.squam.private, x=0.7, y=0.95, gp=gpar(fontsize=15))

## Adeno-private
upset(data=dta[adeno.private, tn$histology == 'CRPC-Ad'], 
      order.by='freq', 
      nsets=sum(tn$histology == 'CRPC-Ad'),
      main.bar.color='#52BCA3',
      sets.bar.color='#52BCA3',
      text.scale=1.9)
grid.text(pct.adeno.private, x=0.7, y=0.95, gp=gpar(fontsize=15))

dev.off()
message(opt$out_file)
