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
## Plot a violin plot of junction burden
libs = c('optparse', 'ggplot2', 'ggbeeswarm', 'ggpubr', 'patchwork')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),       type='character', help="Input junction summary by type"),
  make_option(c("-t", "--tn_file"),       type='character', help="Table specifying which tumor/normal pairs to plot"),
  make_option(c("-m", "--metadata"),      type='character', help="Metadata file"),
  make_option(c("-o", "--out_file"),      type='character', help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor-normal pairs
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t')

if (ncol(tn) == 2) {
  colnames(tn) = c('tumor', 'normal')
} else {
  colnames(tn) = c('tumor','normal','gender')
}
tn$pair_name = paste0(tn$tumor,'--',tn$normal)


## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta$pair_name = paste0(mta$tumor,'--',mta$normal)
tn[, c('site', 'histology')] = mta[match(tn$pair_name, mta$pair_name), c('simplified_site', 'clinical_classification')]


## Read junction burden
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F)
tn$junction_burden = dta$njunc[match(tn$pair_name, dta$sample)]


## Re-bin site
tn$localization = factor(ifelse(tn$site == 'Prostate', 'Prostate', 'Metastatic'), levels=c('Prostate', 'Metastatic'))
tn$localization_liver = ifelse(tn$site == 'Liver', 'Liver', 'non-Liver')



##########
## Plot ##
##########

svg(opt$out_file, width=11, height=6)

hist.col = c(`CRPC-SCC`='#CC61B0', `CRPC-Ad`='#52BCA3')
loc.col = c(`Prostate`='#87C55F', `Metastatic`='#5e92a9')
loc.lvr.col = c(`Prostate`='#87C55F', `Metastatic`='#5e92a9', `Liver`='#F89C74')
loc.lvr.point.col = c(`non-Liver`='#000000', `Liver`='#F89C74')

## Split by histology
p1 = ggplot(tn, aes(x=histology, y=junction_burden, fill=histology)) +
  geom_violin() +
  geom_quasirandom(size=3) +
  scale_fill_manual(values=hist.col) +
  ylab('Junction burden') +
  theme_classic() +
  theme(text=element_text(size=25),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank())

## Split by prostate/met
p2 = ggplot(tn, aes(x=localization, y=junction_burden, fill=localization, col=localization_liver)) +
  geom_violin(col='black') +
  scale_fill_manual(values=loc.col, guide='none') +
  geom_quasirandom(size=3) +
  scale_color_manual(values=loc.lvr.point.col) +
  ylab('Junction burden') +
  theme_classic() +
  theme(text=element_text(size=25),
        legend.position='bottom',
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank())

## Point legend
plt.legend = ggplot(tn, aes(x=localization, y=junction_burden, col=localization_liver)) + 
              geom_point() + 
              scale_color_manual(values=loc.lvr.point.col)
plt.legend = as_ggplot(get_legend(plt.legend))

## Draw
p1 | p2 | plt.legend

dev.off()
message(opt$out_file)
