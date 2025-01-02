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
## Postprocess PhyClone output
libs = c('optparse', 'ggtree', 'treeio', 'data.table', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200)

## https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
large.palette = c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", 
  "#FF7F00", "gold1", "skyblue2", "#FB9A99", 
  "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



## Get arguments
option_list = list(
  make_option(c("-t", "--in_file_tsv"),       type='character', help="TSV output from Phyclone"),
  make_option(c("-p", "--in_file_nwk"),       type='character', help="Newick output from Phyclone"),
  make_option(c("-O", "--out_file_txt"),      type='character', help="Output TSV"),
  make_option(c("-o", "--out_file_pdf"),      type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))


## Read mutations 
mut = fread(opt$in_file_tsv, h=T, stringsAsFactors=F, sep='\t', data.table=F)

## Extract clone CCF and mutation counts
ccf = as.data.frame(tapply(mut$ccf, list(mut$clone_id, mut$sample_id), unique))
ccf$n = as.data.frame(tapply(mut$ccf, list(mut$clone_id, mut$sample_id), length))[, 1, drop=T]
ccf$label = rownames(ccf)
ccf



## Read tree
nwk = read.newick(opt$in_file_nwk)
df = as.data.frame(as_tibble(nwk))


## Assign colors to clones
if (nrow(df) > (length(large.palette) + 1)) {
  message('Too many clones to assign colors')
  node.col = rep('#4989C7', nrow(df))
  names(node.col) = df$label
  node.col['root'] = '#FFFFFF'

} else {
  node.col = rep('#FFFFFF', nrow(df))
  node.col[which(df$label != 'root')] = large.palette[1:(nrow(df)-1)]
  names(node.col) = df$label
  node.col['root'] = '#FFFFFF'

}


## Write data frame with parent-child assignments, ccfs, and colors
df$color = node.col[df$label]
df = merge(x=df, y=ccf, by='label', all.x=T)
write.table(df, opt$out_file_txt, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file_txt)


## Plot tree structure 
pdf(opt$out_file_pdf)

ggplot(nwk, aes(x, y, fill=label)) + 
  geom_tree(layout='roundrect') + 
  geom_nodepoint(shape=21, col='black', size=13, alpha=1) +
  geom_tippoint(shape=21, col='black', size=13, alpha=1) +
  scale_fill_manual(values=node.col) +
  geom_tiplab(hjust=0.8) +
  geom_nodelab(hjust=0.8) +
  theme_tree(legend.position='none')

dev.off()
message(opt$out_file_pdf)
