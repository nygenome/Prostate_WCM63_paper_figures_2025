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
## Plot jabba footprints, including junction support results as annotation track
libs = c('optparse', 'reshape2', 'GenomicRanges', 'gUtils', 'ComplexHeatmap')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Convert a footprint record into GRanges
fp2gr = function(x, id) {
  
  x = gsub('\\+(?=,)|\\-(?=,)','',x, perl=T)
  x = gsub('\\+(?=;)|\\-(?=;)','',x, perl=T)
  x = gsub('\\+$|\\-$','',x, perl=T)
  x = strsplit(x,';|,')
  
  res = list()
  
  for (i in 1:length(x)) {
    
    seqnames = gsub(':.*','',x[[i]])
    start = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[1]))
    end = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[2]))
    
    res[[i]] = GRanges(seqnames, IRanges(start=start, end=end))
  }
  
  res = do.call(c, res)
  res$id = id

  return(res) 
  
}



## Read junction support matrix into gRanges
read.junction.support = function(f, sample_ids) {


  ## Read table
  x = read.csv(opt$junction_support_mtx, h=T, stringsAsFactors=F, check.names=F, row.names=1, sep='\t')
  x = x[, sample_ids]

  ## Split junctions into breakends
  res = strsplit(rownames(x), '|', fixed=T)

  ## For each pair of breakends, expand to GRanges
  res = lapply(res, function(y) makeGRangesFromDataFrame(colsplit(gsub('\\+$|\\-$', '', y), ':|\\-', names=c('chr','start','end'))))

  ## Fill in read counts and ID
  for (i in 1:nrow(x)) {

    mcols(res[[i]]) = x[i, ]

  }

  res = Reduce(c, res)

}



## Compute edit distance between two character vectors
## of equal length
hamming = function(x, y) {

  idx.keep = which(!is.na(x) & !is.na(y))

  res = sum( (x[idx.keep] == 'None') != (y[idx.keep] == 'None') ) / length(x[idx.keep])

  if (is.na(res)) {

    res = 1

  }

  return(res)

}



build.colormap = function(var, col) {
  
  var.col = structure(col, names=var)
  var.col = na.omit(var.col[!duplicated(var)])
  return(var.col)
  
}



events.long = c('Translocation','Inversion', 'Deletion', 'Duplication', 'Chromoplexy', 'Rigma', 'Pyrgo', 'Chromothripsis', 'TIC', 'QRP', 'BFB', 'Double Minute', 'Complex DM', 'Tyfonas', 'Unclassified', 'None')
events.short = c('tra', 'inv', 'del', 'dup', 'chromoplexy', 'rigma', 'pyrgo', 'chromothripsis', 'tic', 'qrp', 'bfb', 'dm', 'cpxdm','tyfonas', 'unclassified_event', 'none')
event.colors  = c(`Translocation`='#B2DF8A', `Inversion`='#c6088d', `Deletion`='#A6CEE3', `Duplication`='#FB9A99',
                  `Chromoplexy`='#33A02C', `Rigma`='#1F78B4', `Pyrgo`='#E31A1C', `Chromothripsis`='darkgreen',
                  `TIC`='orchid4', `QRP`='gold2', `BFB`='firebrick', `Tyfonas`='gray25', `Double Minute`='darkorange',
                  `Complex DM`='darkorange3', `Unclassified`='turquoise',  `None`='#ececec')

event.name.map = structure(events.long, names=events.short)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir"),               type='character', help="Directory holding event footprints"),
  make_option(c("-j", "--junction_support_mtx"), type='character', help="Output of patient-junction-support-by-bed.r"),
  make_option(c("-t", "--tn_file"),              type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-I", "--id_map"),               type='character', help="Map betwen internal and manuscript sample IDs"),
  make_option(c("-p", "--patient"),              type='character', help="Metadata file"),
  make_option(c("-n", "--newick"),               type='character', help="Phylogenetic tree in newick format"),
  make_option(c("-m", "--metadata"),             type='character', help="Sample metadata"),
  make_option(c("-c", "--chr_len"),              type='character', help="Chromosome lengths"),
  make_option(c("-b", "--binsize"),              type='numeric',   help="Size of bins to use. Large == more coarse grained", default=1E4),
  make_option(c("-o", "--out_file"),             type='character', help="Figure output file"))
opt = parse_args(OptionParser(option_list=option_list))



## Read chr lengths
chr.len = read.csv(opt$chr_len, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','end'))
chr.len = chr.len[chr.len$chr %in% paste0('chr', c(1:22, 'X')), ]
chr.len$start = 1
chr.len = makeGRangesFromDataFrame(chr.len)


## Read TN pairs
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, col.names=c('tumor','normal','gender'), sep='\t')
tn$pair_name = paste0(tn$tumor,'--',tn$normal)


if (!is.null(opt$id_map)) { 

  id.map = read.csv(opt$id_map, h=T, stringsAsFactors=F, sep='\t')
  print(head(id.map))
  tn$name = id.map$New.sample.name[match(gsub('-','_',tn$tumor), id.map$Associated.sample)]

} else {

  tn$name = gsub('-(Case-WGS|WGS)','', tn$tumor)

}


## Read metadata and align with tn pairs
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta = mta[match(tn$tumor, mta$tumor), ]
mta$patient[mta$patient_color == '#FFFFFF'] = 'Singleton'

if (!is.null(opt$patient) && !is.null(opt$newick)) {

  dend = labels(phylogram::read.dendrogram(opt$newick))
  dend = dend[dend != 'diploid']
  dend.tumor = gsub('--.*', '', dend)

  tn = tn[tn$patient == opt$patient, ]
  tn = tn[match(dend.tumor, tn$tumor), ]
  tn$sample = dend

}


## Read junction support matrix
js = read.junction.support(f=opt$junction_support_mtx, sample_ids=tn$tumor)
js.mta = granges(js)

## Identify squam- and adeno-private junctions
js.mta$squam.private = as.numeric(apply(mcols(js)[, mta$clinical_classification == 'CRPC-Ad'], 1, function(x) all(x == 0))) ## Private to squam, i.e. no read support in all adenos
js.mta$adeno.private = as.numeric(apply(mcols(js)[, mta$clinical_classification == 'CRPC-SCC'], 1, function(x) all(x == 0))) ## Private to adeno, i.e. no read support in all squams
js.mta$clonal = as.numeric(apply(mcols(js), 1, function(x) all(x > 0))) ## Read support in all 

## Aggregate across tiles
js.tiles = granges(gr.tile(chr.len, opt$binsize))
js.tiles = gr.val(js.tiles, js.mta, val=c('squam.private','adeno.private','clonal'), weighted=F, FUN=sum, default.val=0)

for (i in 1:ncol(mcols(js.tiles))) {
  mcols(js.tiles)[is.na(mcols(js.tiles)[, i]), i] = 0
}


## Read footprints
tn$in_file = paste0(opt$in_dir, '/',tn$pair_name,'/jabba.events.footprints.txt')


fp = lapply(tn$in_file, read.csv, h=T, stringsAsFactors=F, sep='\t')


for (i in 1:length(fp)) {
  
  if (nrow(fp[[i]]) > 0) {
    fp[[i]] = mapply(fp2gr, fp[[i]]$footprint, fp[[i]]$type)
    fp[[i]] = Reduce(c, fp[[i]])
  } else {
    fp[[i]] = GRanges(, id=0)
  }

}

names(fp) = tn$name



## Get events for each sample in opt$binsize bins
tiles = granges(gr.tile(chr.len, opt$binsize))

for (i in 1:length(fp)) {

  fp.i = gr.val(tiles, fp[[i]], 'id')$id
  fp.i[fp.i == ''] = 'none'
  fp.i = event.name.map[gsub(',.*', '', fp.i)]

  mcols(tiles)[, names(fp)[i]] = fp.i

}

seqlevels(tiles) = gsub('^chr', '', seqlevels(tiles))
chr = as.factor(seqnames(tiles))
res = t(as.matrix(mcols(tiles)))



########## 
## Plot ##
##########

pdf(opt$out_file, width=15, height=6)

quantiles = seq(0,1,by=0.2)
tmp = c('grey95', '#c50404')
burden.col = circlize::colorRamp2(c(0, 5), tmp)


## Left annotation
left.ano = ComplexHeatmap::rowAnnotation(`Site`=mta$simplified_site,
                                         `Histology`=mta$clinical_classification,
                                          col=list(`Patient`=build.colormap(var=mta$patient, col=mta$patient_color),
                                                    `Site`=build.colormap(var=mta$simplified_site, col=mta$simplified_site_color),
                                                    `Histology`=build.colormap(var=mta$clinical_classification, col=mta$clinical_classification_color)),
                                          gap=unit(c(0, 0, 0), "mm"),
                                          gp = gpar(col="black"))


## Junction support tracks
ht.clonal.avg = Heatmap(matrix(js.tiles$clonal, nrow=1, dimnames=list('Clonal')),
                    col=burden.col,
                    row_names_side='left',
                    column_split=chr,
                    column_gap = unit(3, "mm"),
                    cluster_rows=F,
                    cluster_columns=F,
                    border=T,
                    heatmap_legend_param = list(title = "Junction burden"))


ht.adeno.avg = Heatmap(matrix(js.tiles$adeno.private, nrow=1, dimnames=list('Adeno-private')),
                    col=burden.col,
                    row_names_side='left',
                    column_split=chr,
                    column_gap = unit(3, "mm"),
                    cluster_rows=F,
                    cluster_columns=F,
                    border=T,
                    show_heatmap_legend=F)


ht.squam.avg = Heatmap(matrix(js.tiles$squam.private, nrow=1, dimnames=list('Squamous-private')),
                    col=burden.col,
                    row_names_side='left',
                    column_split=chr,
                    column_gap = unit(3, "mm"),
                    cluster_rows=F,
                    cluster_columns=F,
                    border=T,
                    show_heatmap_legend=F)


## Main heatmap
ht.main = Heatmap(res,
                col=event.colors,
                column_split=chr,
                row_split=mta$clinical_classification,
                row_names_side='left',
                left_annotation=left.ano,
                column_gap = unit(3, "mm"),
                heatmap_legend_param = list(title = "Event"),
                row_names_gp = gpar(fontsize=8),
                border_gp = gpar(lwd=1),
                border=T,
                cluster_columns=F,
                cluster_rows=T,
                clustering_distance_rows=hamming,
                show_row_names=T)
        

## Draw 
ht.clonal.avg %v% ht.adeno.avg %v% ht.squam.avg %v% ht.main

dev.off()
message(opt$out_file)
