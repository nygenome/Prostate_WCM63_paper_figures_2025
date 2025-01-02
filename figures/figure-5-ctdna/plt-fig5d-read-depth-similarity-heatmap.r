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
## Plot heatmap of cosine similarities between read depth profiles
libs = c('optparse', 'ComplexHeatmap', 'circlize', 'viridisLite')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),   type='character', help="Output of init-compare-read-depth-profiles.r"),
  make_option(c("-m", "--id_map"),    type='character', help="Map between sample IDs"),
  make_option(c("-o", "--out_file"),  type='character', help="Output SVG"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data
dta = read.csv(opt$in_file, row.names=1, h=T, stringsAsFactors=F, sep='\t', check.names=F)
id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
id.map = id.map[id.map$Event_category %in% c('Biopsy', 'cfDNA'), ]
id.map$Associated.sample = gsub('_', '-', id.map$Associated.sample)

rownames(dta) = id.map$New.sample.name[match(rownames(dta), id.map$Associated.sample)]
colnames(dta) = id.map$New.sample.name[match(colnames(dta), id.map$Associated.sample)]


histology = paste0('CRPC-', id.map$Histology[match(rownames(dta), id.map$New.sample.name)])
dta = as.matrix(dta)


## Compare histologies at each timepoint 
cmpr.results = rep(NA, ncol(dta))
for (i in 1:ncol(dta)) {

  cmpr = t.test(x=dta[histology=='CRPC-SCC', i],
                y=dta[histology=='CRPC-Ad', i])

  if (cmpr$p.value < 0.05) {
    
    if (cmpr$estimate['mean of x'] > cmpr$estimate['mean of y']) {
      cmpr.results[i] = 'CRPC-SCC'
    } else {
      cmpr.results[i] = 'CRPC-Ad'
    }


  } else {

    cmpr.results[i] = 'N.S.'

  }

  message(colnames(dta)[i],': ', cmpr$p.value)
  print(cmpr$estimate)

}



## Compare within histologies at timepoints 2/3 
t.test(x=dta[histology=='CRPC-Ad', colnames(dta) == 'WCM63_2'],
       y=dta[histology=='CRPC-Ad', colnames(dta) == 'WCM63_3'])

t.test(x=dta[histology=='CRPC-SCC', colnames(dta) == 'WCM63_2'],
       y=dta[histology=='CRPC-SCC', colnames(dta) == 'WCM63_3'])



## Compare within histologies at timepoints 3/4 
message('\n\n\n')
t.test(x=dta[histology=='CRPC-SCC', colnames(dta) == 'WCM63_3'],
       y=dta[histology=='CRPC-SCC', colnames(dta) == 'WCM63_4'])

t.test(x=dta[histology=='CRPC-Ad', colnames(dta) == 'WCM63_3'],
       y=dta[histology=='CRPC-Ad', colnames(dta) == 'WCM63_4'])



##########
## Plot ##
##########

svg(opt$out_file, width=7, height=8)

nbreaks = 10
col.fun = colorRamp2(quantile(unlist(dta), seq(0, 1, length.out=nbreaks)), plasma(nbreaks))

col.ano = c(`CRPC-Ad`='#52BCA3', `CRPC-SCC`='#CC61B0', `N.S.`='#FFFFFF')
cmpr.results = factor(cmpr.results, levels=names(col.ano))
txt.cutoff = 0.2


## Top annotation 
top.ano = ComplexHeatmap::HeatmapAnnotation(`Closest histology`=cmpr.results,
                                             col=list(`Closest histology`=col.ano),
                                             gp = gpar(col='black'))


Heatmap(dta, 
        name='Cosine similarity',
        col=col.fun, 
        top_annotation=top.ano, 
        row_split=histology, 
        show_row_dend=FALSE,
        cluster_columns=F, 
         heatmap_legend_param = list(ncol=1),
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", dta[i, j]), x, y, gp = gpar(fontsize = 10, fontface='bold', col=ifelse(dta[i, j] < txt.cutoff, 'white', 'black')))}
)

dev.off()
message(opt$out_file)
