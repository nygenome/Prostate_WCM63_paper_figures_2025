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
## Plot jabba track and/or coverage track and/or allele track by a bed file
libs = c('optparse','GenomicRanges','gGnome','gTrack','gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))



## Get arguments
option_list = list(
  make_option(c("-s", "--samples"),         type='character', help="Sample IDs (comma delimited)"),
  make_option(c("-i", "--id_map"),          type='character', help="Map to manuscript IDs"),
  make_option(c("-j", "--jabba_gg"),        type='character', help="Path to jabba.events.rds file"),
  make_option(c("-c", "--coverage"),        type='character', help="Coverage GRanges object (RDS) (optional)"),
  make_option(c("-e", "--metadata"),        type='character', help="Sample metadata (optional)"),
  make_option(c("-f", "--field"),           type='character', help="Field in --coverage file to use when plotting coverage", default='foreground'),
  make_option(c("-b", "--bed"),             type='character', help="BED file with regions for plotting"),
  make_option(c("-p", "--padding"),         type='numeric',   help="Padding for taking subgraph associated with BED intervals"),
  make_option(c("-t", "--junc_table"),      type='character', help="Junction support table with"),
  make_option(c("-m", "--max_junctions"),   type='numeric',   help="Maximum number of junctions to show colors/labels for", default=25),
  make_option(c("-o", "--out_file"),        type='character', help="Figure output file (SVG)"),
  make_option(c("-x", "--out_file_fp"),     type='character', help="Optional output of subgraph footprint"))
opt = parse_args(OptionParser(option_list=option_list))



## Expand arguments
opt$samples = unlist(strsplit(opt$samples, ',',fixed=T))
opt$jabba_gg = unlist(strsplit(opt$jabba_gg, ',',fixed=T))


## Init gene track 
GENCODE = '/gpfs/commons/projects/nepc/analysis/jabba-latest/annotations/gencode.composite.collapsed.rds'

gencode = track.gencode(gencode=GENCODE, 
                          cached.dir=dirname(GENCODE), 
                          cached.path=GENCODE, 
                          cex.label=0.5, 
                          xaxis.cex.label=1, 
                          xaxis.unit=1e6, 
                          xaxis.suffix='MB')
gencode$height=2


if (!is.null(opt$metadata)) {
  mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
  clinical.class = mta$clinical_classification[match(opt$samples, mta$tumor)]
  clinical.class = factor(clinical.class, levels=c('CRPC-Ad', 'CRPC-SCC'))

}


## If we have an ID map, update sample names
if (!is.null(opt$id_map)) {

    id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
    opt$samples = id.map$New.sample.name[match(gsub('-','_',opt$samples), id.map$Associated.sample)]

    if (!is.null(opt$metadata)) {

      ord = order(clinical.class, opt$samples)

      opt$samples = opt$samples[ord]
      opt$jabba_gg = opt$jabba_gg[ord]

    }

}


## Read support table
if (!is.null(opt$junc_table)) {

  junc.support = read.csv(opt$junc_table, h=T, stringsAsFactors=F, row.names=1, sep='\t')
  if ('color' %in% colnames(junc.support)) { 
    junc.support =junc.support[, 'color', drop=F]  
  } else {
    junc.support = NULL
  }

} else {

  junc.support = NULL

}


## Read BED, format intervals
## This is the first mandatory track
bed = rtracklayer::import(opt$bed)


## Add genome graph with edge colors aligned to junction support heatmap
jabba.gt = vector(mode='list', length=length(opt$samples))
jabba.fp = vector(mode='list', length=length(opt$samples))
jabba.ag.gt = vector(mode='list', length=length(opt$samples))

for (i in 1:length(opt$samples)) {
  
  ## Read genome graph
  jabba.gg = readRDS(opt$jabba_gg[i])
  jabba.gg$nodes$mark(col='gray')
  
  ## Order gTracks to match junction support table 
  ## Mark edges by junction ID, don't bother if we have too many junctions
  ## as many colors will make the plot impossible to decipher
  if (!is.null(junc.support) && nrow(junc.support) < opt$max_junctions) {

    ## Reset node/edge coloring
    jabba.gg$nodes$mark(col='gray')
    jabba.gg$edges[type=='ALT']$mark(col='gray50')

    ## Get string representation of junctions
    junc.gr = grl.unlist(jabba.gg$junctions[type=='ALT']$grl)
    junc.str = gr.string(junc.gr)
    junc.str = sapply(1:max(junc.gr$grl.ix), function(i) paste(junc.str[junc.gr$grl.ix==i], collapse='|'))

    for (j in 1:nrow(junc.support)) {
      
      junc.id = rownames(junc.support)[j]
      junc.id.rev = paste(rev(unlist(strsplit(junc.id,'\\|'))), collapse='|')
      junc.col = junc.support$color[j]
      
      ## Determine which edge to mark
      edge.id.sel = unique(junc.gr$edge.id)[junc.str == junc.id | junc.str == junc.id.rev]
      
      if (length(edge.id.sel) == 0) {
        next
      }
      
      jabba.gg$edges[edge.id == eval(edge.id.sel)]$mark(col=junc.col)
      
    }

  }
  
  ## Build gtrack
  jabba.fp[[i]] = jabba.gg$copy$subgraph(bed, opt$padding)$footprint
  jabba.gt[[i]] = jabba.gg$gtrack(name=opt$samples[i])

}


## Take union of footprints for plotting region
jabba.fp = do.call(gr.reduce, jabba.fp)


## Combine gTracks
for (i in length(opt$samples):1) {
  
  if (i == length(opt$samples)) {
    plt.gtrack = jabba.gt[[i]]
  } else {
    plt.gtrack = c(plt.gtrack, jabba.gt[[i]])
  }
  
}



##########
## Plot ## 
##########

## Scale padding based on the largest window 
plt.windows = reduce(c(jabba.fp,bed))
plt.windows = plt.windows + (max(width(plt.windows)) * 0.1)

## Scale plot area and y-axis labels based on the number of samples we're plotting 
height.scale.factor = 2 * (1 - (length(opt$samples) / 4) )
cex.y.scale.factor = 1.25 - (0.0625 * length(opt$samples))

## Scale x-axis labels based on the number of windows we're plotting 
cex.x.scale.factor = 1.25 - (0.01 * length(plt.windows))

svg(opt$out_file, width=10, height=10 - height.scale.factor)

plot(c(gencode, plt.gtrack), 
     plt.windows, 
     cex.ylabel=cex.y.scale.factor, 
     yaxis.cex=cex.y.scale.factor, 
     cex.xlabel=cex.x.scale.factor, 
     xaxis.cex=cex.x.scale.factor,
     xaxis.cex.tick=cex.x.scale.factor)
  
dev.off()

message(opt$out_file)
