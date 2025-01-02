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
## Check junction support across a cohort specific loci, 
## or across the entire genome if no loci given
libs = c('optparse', 'GenomicRanges', 'gGnome', 'gUtils', 'skitools', 'RSeqLib', 'Biostrings', 'RColorBrewer', 'rtracklayer')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Read a bed file into granges object
read.bed = function(f) {
  
  x = read.csv(f, h=F, stringsAsFactors=F, sep='\t')

  if (ncol(x) == 3) {
    colnames(x) = c('chr','start','end')
  } else {
    colnames(x) = c('chr','start','end', 'id')
  }

  x = makeGRangesFromDataFrame(x, keep.extra=T)
  
  return(x)
  
}



loadJunctions = function(f, gr) {
  
  gg = gG(jabba=f)
  
  ## Trim to region of interest
  junc = gg$junctions[type == 'ALT'] 
  
  if (!is.null(gr)) {
    junc = junc[sapply(junc$grl, function(x) any(x %^% gr))]  
  }
  
  print(junc)
  
  return(junc)
  
}



## Get arguments
option_list = list(
  make_option(c("-j", "--jba_dir"),            type='character', help="Directory holding JaBbA results"),
  make_option(c("-a", "--bam_dir"),            type='character', help="Top level project dir with BAMs i.e. bam_dir/Sample_*/analysis/*final.bam"),
  make_option(c("-t", "--tn_file"),            type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-b", "--bed"),                type='character', help="Comma-delimited list of bed files to select junctions"),
  make_option(c("-p", "--pad"),                type='numeric',   help="Optional padding for BED intervals (default=0)", default=0),
  make_option(c("-s", "--samples"),            type='character', help="Comma-delimited list of tumor IDs (will automatically pull normals as controls)"),
  make_option(c("-d", "--additional_samples"), type='character', help="Additional samples that we have BAMs for, but might not be part of the TN list"),
  make_option(c("-g", "--genome"),             type='character', help="Reference genome FASTA with BWA index in the same dir"), 
  make_option(c("-o", "--out_file"),           type='character', help="Output file"))
opt = parse_args(OptionParser(option_list=option_list))



JUNCTION_PADDING = 2000
REF_PADDING = 1E4
MERGE_PADDING = 300


## Read BWA index, reference genome as DNAStringSet
bwa = RSeqLib::BWA(opt$genome)
ref = readDNAStringSet(opt$genome)


## Read optional BED(s)
if (!is.null(opt$bed)) {
  
  bed = lapply(unlist(strsplit(opt$bed,',')), read.bed)
  bed = gr.reduce(Reduce(f=c, bed)) + opt$pad
  
} else {
  
  bed = NULL
  
}


## Optionally adjust footprints for simple events that span a wide footprint 
## But are only composed of 1-2 junctions
## Otherwise any junctions inside come along for the ride 
if (!is.null(opt$bed) && 'id' %in% colnames(mcols(bed))) {

  idx.simple = grep('(dup|del|inv)', bed$id)
  if (length(idx.simple) > 0) {
    bed = c(gr.start(bed), gr.end(bed)) + 1E3
  }

}


## Read tumor normal pairing file
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)
samples = unique(c(tn$tumor, tn$normal))  

jba.files = paste0(opt$jba_dir,'/',tn$pair_name,'/jabba.simple.gg.rds')
bam.files = paste0(opt$bam_dir,'/Sample_',samples,'/analysis/',samples,'.final.bam')


## Take union of junctions with no padding
juncs = lapply(jba.files, loadJunctions, gr=bed)
junc.count = sapply(juncs, length)
juncs$pad = MERGE_PADDING
junc.union = do.call(gGnome::merge, juncs[junc.count > 0])


## Take reads around junctions of interest
bam.window = gUtils::streduce(junc.union$breakpoints, pad=JUNCTION_PADDING)

## Create dataframe to hold result
juncs.gr = grl.unlist(junc.union$grl)
juncs.str = gr.string(juncs.gr)
juncs.str = sapply(1:max(juncs.gr$grl.ix), function(i) paste(juncs.str[juncs.gr$grl.ix==i], collapse=('|')))

res = as.data.frame(matrix(0, nrow=length(juncs.str), ncol=length(samples)))
rownames(res) = juncs.str
colnames(res) = samples

bin.mtx = as.data.frame(junc.union$dt)



## If one sample has no junctions in this region, it won't show up here
## Fill in blank columns
if (any(junc.count == 0)) {
  
  missing = tn$tumor[which(junc.count == 0)]
  present = tn$tumor[which(junc.count > 0)]
  

  bin.mtx = bin.mtx[, grep('^seen\\.by\\.ra', colnames(bin.mtx)) , drop=F]
  colnames(bin.mtx) = paste0('seen.by.', present)
  
  bin.mtx[, paste0('seen.by.', missing)] = rep(FALSE, nrow(bin.mtx))
  rownames(bin.mtx) = rownames(res)
  
} else {
  
  bin.mtx = bin.mtx[, paste0('seen.by.ra',1:nrow(tn)), ]
  colnames(bin.mtx) = paste0('seen.by.', tn$tumor)
  rownames(bin.mtx) = rownames(res)
    
}



## Examine BAMs for read support
for (i in 1:length(bam.files)) {


  bam = bamUtils::read.bam(bam.files[i], intervals=bam.window, pairs.grl=F, tag="AS")
  message('Read BAM for ', bam.files[i])

  ## The seqlengths that come out of read.bam can be inconsistent depending on the interval chosen
  ## So recover the correct ones from the header
  seqlength.update = seqlengths(Rsamtools::BamFile(bam.files[i]))
  seqlengths(bam) = seqlength.update[names(seqlengths(bam))]

  ## Compute read support for each junction
  res.i = skitools::junction.support(reads=bam,
                                     junctions=junc.union,
                                     bwa=bwa,
                                     ref=ref,
                                     realign=F,
                                     pad=JUNCTION_PADDING,
                                     pad.ref=REF_PADDING)

  read.support = table(res.i$junction.id)
  res[as.numeric(names(read.support)), i] = unname(read.support)
  
}

## Combine read support and presence/absence info
res = cbind(res, bin.mtx)



## Assign a unique color to each junction
## Right now only works with up to 17 junctions, any more would probably be 
## too busy anyways
cols = c(brewer.pal(8, 'Dark2')[-c(4,5,7,8)], brewer.pal(9, 'Set1'), brewer.pal(9, 'Pastel1'))
if (nrow(res) < length(cols)) {
  res$color = cols[1:nrow(res)]
} else {
  res$color = sample(cols, nrow(res), replace=T)
}

  
## Write result
write.table(res, opt$out_file, row.names=T, col.names=T, quote=F, sep='\t')
message(opt$out_file)
