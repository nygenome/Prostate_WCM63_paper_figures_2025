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
## Plot a junction burden heatmap
libs = c('optparse', 'ComplexHeatmap', 'circlize', 'dendextend', 'ggplot2', 'RColorBrewer')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



build.colormap = function(var, col) {
  
  var.col = structure(col, names=var)
  var.col = na.omit(var.col[!duplicated(var)])
  return(var.col)
  
}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),       type='character', help="Input junction summary by type"),
  make_option(c("-t", "--tn_file"),       type='character', help="Table specifying which tumor/normal pairs to plot"),
  make_option(c("-m", "--metadata"),      type='character', help="Metadata file"),
  make_option(c("-M", "--id_map"),        type='character', help="Map to manuscript IDs (optional)"),
  make_option(c("-k", "--k"),             type='numeric',   help="Number of clusters to select", default=6),
  make_option(c("-o", "--out_file"),      type='character', help="Output PDF"))
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
full.sample.count = nrow(mta)
mta$pair_name = paste0(mta$tumor,'--',mta$normal)
mta = mta[mta$pair_name %in% tn$pair_name, ]

## Read data, select tumor-normal pairs, convert to matrix
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F)


## Align data/metadata
rownames(dta) = dta$sample
rownames(mta) = mta$pair_name
pairs.keep = intersect(dta$sample, mta$pair_name)
mta = mta[pairs.keep, ]
dta = dta[pairs.keep, setdiff(colnames(dta), 'sample')]

njunc = dta$njunc
dta$dm = dta$dm + dta$cpxdm 
dta = as.matrix(dta[, setdiff(colnames(dta), c('unclassified', 'njunc', 'cpxdm'))])


## Take ln(x+1) and transpose
dta = t(log1p(dta))


# Separate simple/complex events
simple.events = c('del', 'dup', 'inv', 'tra', 'invdup')
simple.complex = ifelse(rownames(dta) %in% simple.events, 'Simple', 'Complex')


## Update event names
rownames(dta) = gsub('invdup','Inverted Duplication',rownames(dta))
rownames(dta) = gsub('tra','Translocation',rownames(dta))
rownames(dta) = gsub('inv','Inversion',rownames(dta))
rownames(dta) = gsub('del','Deletion',rownames(dta))
rownames(dta) = gsub('dup','Duplication',rownames(dta))
rownames(dta) = gsub('chromoplexy','Chromoplexy',rownames(dta))
rownames(dta) = gsub('rigma','Rigma',rownames(dta))
rownames(dta) = gsub('pyrgo','Pyrgo',rownames(dta))
rownames(dta) = gsub('chromothripsis','Chromothripsis',rownames(dta))
rownames(dta) = gsub('tic','TIC',rownames(dta))
rownames(dta) = gsub('qrp','QRP',rownames(dta))
rownames(dta) = gsub('cpxdm','Complex DM',rownames(dta))
rownames(dta) = gsub('bfb','BFB',rownames(dta))
rownames(dta) = gsub('dm','Double Minute',rownames(dta))
rownames(dta) = gsub('tyfonas','Tyfonas',rownames(dta))


## Remove normal from sample names
colnames(dta) = gsub('--.*','',colnames(dta))


## Optionally update names
if (!is.null(opt$id_map)) {

    id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
    colnames(dta) = id.map$New.sample.name[match(gsub('-','_',colnames(dta)), id.map$Associated.sample)]

} else {

  colnames(dta) = gsub('-Case-WGS','',colnames(dta))

}



##############
## Plotting ##
##############

mta$patient_name = mta$patient
mta$patient_name[mta$patient_color == '#FFFFFF'] = 'Singleton'

site.col = build.colormap(var=mta$simplified_site, col=mta$simplified_site_color)
patient.col = build.colormap(var=mta$patient_name, col=mta$patient_color)
clinical.col = build.colormap(var=mta$clinical_classification, col=mta$clinical_classification_color)
 
top.ano = ComplexHeatmap::HeatmapAnnotation(`Junction Count`=ComplexHeatmap::anno_barplot(njunc, baseline=0, gp=gpar(fill='black')),
                                            `Site`=mta$simplified_site,
                                            `Histology`=mta$clinical_classification,
                                            col=list(`Patient`=patient.col,
                                                      `HRD`=c(`TRUE`='#1f1f1f', `FALSE`='#e8e8e8'),
                                                      `Site`=site.col,
                                                      `Histology`=clinical.col),
                                            gap=unit(c(2, 0, 0, 0, 0), "mm"),
                                            gp = gpar(col="black"),
                                            show_legend=T)


quantiles = seq(0,1,by=0.2)
tmp = brewer.pal(n = length(quantiles), name = "YlOrRd")
tmp = c('grey95', tmp)
burden.col = circlize::colorRamp2(c(0, quantile(dta[dta > 1], quantiles)), tmp)


## Compute dendrogram
k = opt$k
dend = hclust(dist(t(dta), 'euclidean'), method='complete')
dend = dendextend::color_branches(dend, k=k)


## Main burden heatmap
main.heatmap = ComplexHeatmap::Heatmap(dta,
                                      col=burden.col,
                                      border=T,
                                      top_annotation=top.ano,
                                      row_names_side='left',
                                      column_title=NULL,
                                      heatmap_legend_param = list(title = "ln(burden)"),
                                      show_column_names=T,
                                      show_row_names=T,
                                      cluster_rows=T,
                                      cluster_columns=dend,
                                      column_split=k,
                                      row_split=simple.complex,
                                      show_heatmap_legend=T,
                                      column_names_gp = gpar(fontsize = 14))


## Draw heatmap
width.scaling = min((8 * ncol(dta))/full.sample.count, 1)
pdf(opt$out_file, width=19 * width.scaling, height=10)

draw(main.heatmap)

dev.off()
message(opt$out_file)
