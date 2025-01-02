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
## Compare junction support within events 
libs = c('optparse', 'ComplexHeatmap' , 'circlize')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



option_list = list(
  make_option(c("-i", "--in_file"),   type='character', help="Output of patient-junction-support-by-bed.r"),
  make_option(c("-m", "--metadata"),  type='character', help="Sample metadata"),
  make_option(c("-M", "--id_map"),    type='character', help="Map to manuscript IDs"),
  make_option(c("-o", "--out_file"),  type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



MIN_READ_SUPPORT = 1 
MIN_CONSERVATION = 0.75


## Read data
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, check.names=F, sep='\t')
dta = dta[, grep('^seen', colnames(dta), invert=T, value=T)]
all.samples = setdiff(colnames(dta), c('id', 'junction'))


## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
clinical.class = mta$clinical_classification[match(all.samples, mta$tumor)]



############################
## Summarize conservation ##
############################

## Which junctions have at least one read of support? 
has.support = as.data.frame(apply(dta[, all.samples], 2, `>=`, 1))
has.support = split(has.support, dta$id)

res.cons = data.frame(id=unique(dta$id), num_squam=NA, num_adeno=NA, histology=NA)

for (i in 1:nrow(res.cons)) {
  
  res.cons$num_squam[i] = sum(apply(has.support[[res.cons$id[i]]][, clinical.class=='CRPC-SCC'], 
                              2, 
                              function(x) (sum(x) / length(x)) >= MIN_CONSERVATION))

  res.cons$num_adeno[i] = sum(apply(has.support[[res.cons$id[i]]][, clinical.class=='CRPC-Ad'], 
                              2, 
                              function(x) (sum(x) / length(x)) >= MIN_CONSERVATION))

  if (res.cons$num_squam[i] > 0 & res.cons$num_adeno[i] > 0) {
    res.cons$histology[i] = 'Both'
  } else if (res.cons$num_squam[i] > 0 & res.cons$num_adeno[i] == 0) {
    res.cons$histology[i] = 'CRPC-SCC'
  } else if (res.cons$num_squam[i] == 0 & res.cons$num_adeno[i] > 0) {
    res.cons$histology[i] = 'CRPC-Ad'
  } 

}

## If event isn't conserved in any samples, file under "Both"
res.cons$histology[res.cons$num_squam == 0 & res.cons$num_adeno == 0] = 'Both'

message('In squamous but not in adeno')
res.cons[res.cons$num_squam > 0 & res.cons$num_adeno == 0, ]

message('\n\nIn adeno but not in squamous')
res.cons[res.cons$num_squam == 0 & res.cons$num_adeno > 0, ]


table((res.cons$num_adeno + res.cons$num_squam) == 12)
res.cons$histology[(res.cons$num_adeno + res.cons$num_squam) == 12] = 'All'
dim(res.cons)

table(res.cons$histology)

dta$histology = res.cons$histology[match(dta$id, res.cons$id)]



#####################
## Support heatmap ##
#####################

read.support = as.matrix(dta[, all.samples])

## Color scale 
col.fun = colorRamp2(c(0,0.99, quantile(read.support, 0.95)), c("gray90","white", "red"))
lgd = Legend(labels = 'Not in Jabba model', title = "", type = "points",
             pch='x', legend_gp = gpar(col = '#000000'), background = "white")


## Align data
clinical.cols = c(`CRPC-Ad`='#52BCA3', `CRPC-SCC`='#CC61B0')
clinical.class = factor(clinical.class, c('CRPC-Ad', 'CRPC-SCC'))
ord = order(clinical.class)
 
read.support = read.support[, ord]
clinical.class = clinical.class[ord]


## Top annotation 
top.ano = ComplexHeatmap::HeatmapAnnotation(`Histology`=clinical.class,
                                             col=list(`Histology`=clinical.cols),
                                             gp = gpar(col='black'))



##########
## Plot ##
##########

pdf(opt$out_file, width=4, height=15)

junc.per.event = table(dta$id)
dta$complex = dta$id %in% names(junc.per.event)[junc.per.event > 2]

dta$id = gsub('PM63_union_footprint_', '', dta$id, fixed=T)


for (i in c('All', 'Both', 'CRPC-SCC', 'CRPC-Ad')) {

  plt = Heatmap(read.support[which(dta$histology == i & !dta$complex), ],
        col=col.fun,
        border=T,
        column_split=clinical.class,
        row_split=dta$id[dta$histology == i & !dta$complex], 
        top_annotation=top.ano,
        column_title=paste0(i,' - "Simple"'),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        heatmap_legend_param = list(title = "Read support"),
        show_row_names=F,
        show_row_dend=F,
        cluster_rows=T,
        cluster_columns=F)

  draw(plt, annotation_legend_side='bottom')

  plt = Heatmap(read.support[which(dta$histology == i & dta$complex), ],
        col=col.fun,
        border=T,
        column_split=clinical.class,
        row_split=dta$id[dta$histology == i & dta$complex], 
        top_annotation=top.ano,
        column_title=paste0(i,' - "Complex"'),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        heatmap_legend_param = list(title = "Read support"),
        show_row_names=F,
        show_row_dend=F,
        cluster_rows=T,
        cluster_columns=F)

  draw(plt, annotation_legend_side='bottom')

}


dev.off()
message(opt$out_file)
