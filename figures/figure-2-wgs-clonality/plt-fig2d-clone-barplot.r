#!/nfs/sw/R/R-4.0.0/bin/Rscript
#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2025) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper & Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################
## Plot per-sample clone percentages
libs = c('optparse', 'reshape2', 'ggplot2', 'patchwork')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)
set.seed(1)



## Recursively constrain ccfs so that total CCF of observed clones sums to 1
## For each node,
##   if it's a leaf, return the current value
##   if it's not a leaf, subtract the sum of the children's value from this node 
constrain.ccfs = function(x) {

  ## Get sorted list of parents, remove root as it doesn't have any variants
  parents = setdiff(unique(x$parent), 'root')
  parents = as.character(sort(as.numeric(parents)), decreasing=F)

  ## For each parent, working our way from root to tip, 
  ## subtract child frequencies
  for (p in parents) {
    x$ccf[x$node == p] = x$ccf[x$node == p] - sum(x$ccf[x$parent == p])
  }

  return(x$ccf)

}



## Normalized entropy 
## Adapted from 
## https://mc-stan.org/posterior/reference/entropy.html
## https://github.com/stan-dev/posterior/blob/ed244711a9b558c356aa51bae82e3ab9ab919c88/R/discrete-summaries.R#L54
normalized.entropy <- function(x) {

  if (anyNA(x)) return(NA_real_)

  p = x 
  n = length(p)

  if (n == 1) {
    out = 0
  } else {
    p = p[p > 0]
    out = -sum(p * log(p)) / log(n)
  }

  out

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Output of analysis/02-clonal-deconvolution/init-post-process-phyclone.r. Table of clone relationships and CCFs"),
  make_option(c("-m", "--id_map"),   type='character', help="Map between internal and manuscript IDs"),
  make_option(c("-o", "--out_file"), type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))


## Read ID map
id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
id.map = id.map[id.map$Event_category == 'Biopsy', ]
id.map$Associated.sample = gsub('_1_Case.*|_2_Case_WGS', '', id.map$Associated.sample)


## Read data 
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
dta = dta[dta$label != 'root', !colnames(dta) %in% c('n')]


## Update node names 
dta$parent = dta$label[match(dta$parent, dta$node)]
dta$node = dta$label
dta = dta[, colnames(dta) != 'label']
dta$parent[is.na(dta$parent)] = 'root'


## Constrain displayed CCFs based on tree
sample.ids = setdiff(colnames(dta), c('parent', 'node', 'color'))

for (s in sample.ids) {

  dta.s = dta[, c('parent', 'node', s)]
  colnames(dta.s)[colnames(dta.s) == s] = 'ccf'
  dta[, s] = constrain.ccfs(x=dta.s)

}


## Compute entropy 
dta.entropy = data.frame(Sample=sample.ids, Entropy=NA)
for (i in 1:nrow(dta.entropy)) {
  dta.entropy$Entropy[i] = normalized.entropy(dta[, dta.entropy$Sample[i]])
}


## Reshape for plotting
col = dta$color
names(col) = dta$node
dta$node = factor(dta$node, levels=as.character(sort(as.numeric(dta$node), decreasing=F)))
dta = melt(dta[, c('node', sample.ids)], 
           id.var='node', 
           variable.name='Sample', 
           value.name='CCF')
dta = dta[dta$CCF > 0, ]


## Round CCFs and drop hidden nodes
dta$Clone = droplevels(dta$node)
dta$CCF = round(dta$CCF, 4)


## Rename samples
dta$Sample = id.map$New.sample.name[match(dta$Sample, id.map$Associated.sample)]
dta.entropy$Sample = id.map$New.sample.name[match(dta.entropy$Sample, id.map$Associated.sample)]


## Set ordering manually
message('Manually setting sample ordering')
dta$Sample = factor(gsub('WCM63_', '', dta$Sample), levels=c('A', 'D', 'E', 'C', 'F', 'G', 'B', 'H', 'I', 'L', 'J', 'K'))
dta.entropy$Sample = factor(gsub('WCM63_', '', dta.entropy$Sample), levels=c('A', 'D', 'E', 'C', 'F', 'G', 'B', 'H', 'I', 'L', 'J', 'K'))


## Set histology manually
message('Manually setting histology')
dta$Histology = ifelse(dta$Sample %in% c('A', 'D', 'C', 'E', 'F', 'G'), 'Adeno', 'Squamous')
dta.entropy$Histology = ifelse(dta.entropy$Sample %in% c('A', 'D', 'C', 'E', 'F', 'G'), 'Adeno', 'Squamous')
dta.entropy$Primary = ifelse(dta.entropy$Sample %in% c('A', 'D'), 'Primary', 'Metastatic')


message('Comparison of normalized Shannon entropy between primary/mets: ')
t.test(Entropy ~ Primary, data=dta.entropy, paired=F)



##########
## Plot ##
##########

svg(opt$out_file, width=6, height=4)

plt.main = ggplot(dta, aes(x=Sample, y=CCF, fill=Clone)) +
            geom_bar(stat='identity', position='stack') +
            facet_wrap(. ~ Histology, nrow=1, scale='free_x') +
            scale_fill_manual(values=col) +
            theme_bw() + 
            theme(text=element_text(size=15),
                  strip.background=element_blank(),
                  strip.text=element_blank(), 
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(), 
                  axis.title.x=element_blank(),
                  legend.position='right')

plt.entropy = ggplot(dta.entropy, aes(x=Sample, y=Entropy)) + 
                geom_bar(stat='identity') +
                facet_wrap(. ~ Histology, nrow=1, scale='free_x') +
                theme_bw() + 
                theme(text=element_text(size=15),
                      strip.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(), 
                      axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.x=element_blank(),
                      legend.position='right')

plt.entropy + plt.main + plot_layout(nrow=2, ncol=1, heights=c(2, 8))

dev.off()
message(opt$out_file)
