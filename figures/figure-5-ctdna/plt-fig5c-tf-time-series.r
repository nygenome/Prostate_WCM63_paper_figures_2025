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
## Time series plots of clone-level detection rates in plasma
libs = c('optparse', 'reshape2', 'ggplot2', 'patchwork', 'scales')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file_tf"),   type='character', help="TF estimates"),
  make_option(c("-m", "--id_map"),       type='character', help="ID map with information about day of sampling"),
  make_option(c("-t", "--in_file_tree"), type='character', help="PhyClone tree table"),
  make_option(c("-o", "--out_file"),     type='character', help="Output SVG"))
opt = parse_args(OptionParser(option_list=option_list))


## Read data
tf = read.csv(opt$in_file_tf, h=T, stringsAsFactors=F, sep='\t')
clones = c('all', as.character(sort(as.numeric(setdiff(unique(tf$clone), 'all')))))
tf$clone = factor(tf$clone, levels=clones)

tree = read.csv(opt$in_file_tree, h=T, stringsAsFactors=F, sep='\t')
tree = tree[tree$label != 'root', ]


## Determine which clones are shared across >1 samples 
sample.cols = setdiff(colnames(tree), c('label', 'parent', 'node', 'n', 'color'))
tree$multiple_samples = apply(tree[, sample.cols], 1, function(x) sum(x > 0) > 1)
tree = tree[tree$multiple_samples, ]


## Set 0 to very low value for log scaling
tf$detection_rate.ci95.lower[tf$detection_rate.ci95.lower == 0] = 1E-6


## Exclude patient-wide detection rate 
tf = tf[tf$clone != 'all', ]


## Update clone IDs
tf = tf[tf$clone %in% tree$label, ]
tf$clone = factor(paste0('C', tf$clone), levels=paste0('C', sort(unique(tf$clone))))


## Assign colors
col = tree$color
names(col) = paste0('C', as.character(tree$label))
col = col[as.character(tf$clone)]


## Annotate with day of sampling
id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
tf$day = id.map$Day[match(tf$sample_id, id.map$New.sample.name)]


## Get cutoffs
tf.cutoff = tf[!duplicated(tf$clone), c('clone', 'tf_cutoff')]



##########
## Plot ##
##########

svg(opt$out_file, height=6, width=8)

ggplot(tf, aes(x=day, y=detection_rate, col=clone, group=clone)) +
  geom_line() + 
  geom_point() + 
  geom_ribbon(aes(ymin=detection_rate.ci95.lower, 
                 ymax=detection_rate.ci95.upper), 
              linetype=2, 
              alpha=0.1) +
  facet_wrap(. ~ clone) +
  geom_hline(data=tf.cutoff, aes(yintercept=tf_cutoff), linetype='dashed', alpha=1) +
  scale_color_manual(values=col) +
  scale_x_continuous(name='Day', breaks=sort(unique(tf$day))) +
  scale_y_continuous(name='Plasma detection rate', limits=c(0,max(tf$detection_rate))) +
  theme_bw() +
  theme(text=element_text(size=12),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        strip.background=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7.5),
        legend.position='none')

dev.off()
message(opt$out_file)
