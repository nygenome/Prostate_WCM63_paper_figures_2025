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
## Plot heatmap of read support per junction
libs = c('optparse', 'ComplexHeatmap', 'circlize')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



option_list = list(
  make_option(c("-i", "--in_file"),       type='character', help="Output of patient-junction-support-by-bed.r"),
  make_option(c("-m", "--metadata"),      type='character', help="Sample metadata"),
  make_option(c("-M", "--id_map"),        type='character', help="Map to manuscript IDs"),
  make_option(c("-j", "--max_junctions"), type='numeric',   help="Maximum number of junctions to show colors/labels for", default=25),
  make_option(c("-o", "--out_file"),      type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')

## Read data
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, check.names=F, row.names=1, sep='\t')


## Get colors and assign junction IDs
junc.colors = dta$color
names(junc.colors) = make.unique(rep(LETTERS, ceiling(nrow(dta) / length(LETTERS)))[1:nrow(dta)], sep='')
dta = dta[, colnames(dta) != 'color']
rownames(dta) = names(junc.colors)

all.samples = grep('^seen', colnames(dta), invert=T, value=T)

## Add dummy "seen by" columns for any normal samples
normals = all.samples[!paste0('seen.by.', all.samples) %in% colnames(dta)]
dta[, paste0('seen.by.', normals)] = FALSE

clinical.class = mta$clinical_classification[match(all.samples, mta$tumor)]
clinical.class[is.na(clinical.class)] = 'Normal'



###################
## Build heatmap ##
###################

read.support = as.matrix(dta[, all.samples])

seen.by = dta[, paste0('seen.by.', all.samples)]
seen.by = apply(seen.by, 2, ifelse, '', 'x')

## Don't add jabba info for non-tumor samples
seen.by[, paste0('seen.by.',normals)] = ''

## Split by tumor/non-tumor
col.split = ifelse(colnames(read.support) %in% normals, ' ', 'Tumors')
col.split = factor(col.split, levels=c('Tumors', ' '))


## If we have an ID map, update sample names
if (!is.null(opt$id_map)) {

  id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
  new.col.names = id.map$New.sample.name[match(gsub('-','_',colnames(read.support)), id.map$Associated.sample)]
  new.col.names[is.na(new.col.names)] = 'Normal'
  colnames(read.support) = new.col.names

}


## Color scale 
col.fun = colorRamp2(c(0,0.99, quantile(read.support, 0.95)), c("gray90","white", "red"))
lgd = Legend(labels = 'Not in Jabba model', title = "", type = "points",
             pch='x', legend_gp = gpar(col = '#000000'), background = "white")


## Align data
clinical.cols = c(`CRPC-Ad`='#52BCA3', `CRPC-SCC`='#CC61B0', `Normal`='#FFFFFF')
clinical.class = factor(clinical.class, c('CRPC-Ad', 'CRPC-SCC', 'Normal'))
ord = order(clinical.class, colnames(read.support))

seen.by = seen.by[, ord]
col.split = col.split[ord]
read.support = read.support[, ord]
clinical.class = clinical.class[ord]



##########
## Plot ##
##########

## Left-hand annotation
left.ano = ComplexHeatmap::rowAnnotation(`ID`=names(junc.colors),
                                         col=list(`ID`=junc.colors),
                                         show_legend=FALSE,
                                         show_annotation_name=FALSE,
                                         gp = gpar(col="black"))

## Top annotation 
top.ano = ComplexHeatmap::HeatmapAnnotation(`Histology`=clinical.class,
                                             col=list(`Histology`=clinical.cols),
                                             gp = gpar(col='black'))


## Update some display features depending on how many junctions we have 
if (nrow(read.support) > opt$max_junctions) {

  left.ano = NULL
  cell_fun = NULL

} else {

  cell_fun = function(j, i, x, y, w, h, col) grid.text(seen.by[i, j], x, y)

}



## Handle row grouping
row.split = rep(NA, nrow(read.support))
row.split.read.support.cutoff = 5

for (i in 1:nrow(read.support)) {

  rs.i = unlist(read.support[i, which(colnames(read.support) != 'Normal')])
  cc.i = as.character(clinical.class[which(colnames(read.support) != 'Normal')])

  if (all(rs.i > 0)) {

    row.split[i] = 'Clonal'

  } else if (length(unique(cc.i[rs.i > row.split.read.support.cutoff])) > 1) {
    
    row.split[i] = 'Mixed'

  } else {

    row.split[i] = unique(cc.i[rs.i > row.split.read.support.cutoff])

  }

}

row.split = factor(row.split, levels=c('Clonal', 'CRPC-Ad', 'CRPC-SCC', 'Mixed'))
levels(row.split) = c('I. Clonal', 'II. CRPC-Ad', 'III. CRPC-SCC', 'Mixed')


plt = Heatmap(read.support,
        col=col.fun,
        border=T,
        column_split=clinical.class,
        row_split=row.split,
        row_title_gp = gpar(fontsize = 12),
        left_annotation=left.ano,
        top_annotation=top.ano,
        column_title=NULL,
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        heatmap_legend_param = list(title = "Read support"),
        show_row_names=FALSE,
        row_names_side='left',
        cluster_rows=F,
        cluster_columns=F,
        cell_fun = cell_fun)

pdf(opt$out_file, width=4.2, height=7.5 + 0.05*nrow(read.support))

draw(plt, row_title='Junctions', newpage=F, annotation_legend_list=list(lgd), annotation_legend_side='bottom')

dev.off()
message(opt$out_file)
