library(DNAcopy)
options(scipen=500) # turn off scientific notation; default is 0
# Mon 01 Nov 2021 05:22:31 PM PDT
# script to run DNAcopy
# DNAcopy takes log ratio data (similar to CopyNumber)
# use the corrected and normalized data from HMMcopy
# following steps from https://bioconductor.org/packages/release/bioc/vignettes/DNAcopy/inst/doc/DNAcopy.pdf

# run through bedtools after to get to a windowed file (this script outputs segments)
# ex. bedtools intersect -a simu_cancer_cell_0.hg19_lite.dnaCopy.bed -b /space/s2/sandra/hg19/hg19_lite.250kb_windows | sort -k1,1V -k2,2n
# also need to scale later to get copy number calls (a la CopyNumber)

# Arg 1: path to tumor file csv basename (auto appends .csv to basenames)
# Arg 2: /path/to/output/dir

args <- commandArgs(trailingOnly=T)
tumorBase <- args[1] # cell name
outputBase <- args[2]

dat <- read.table(paste0(tumorBase, ".csv"), sep=",", stringsAsFactors=F, header=T)
sampName <- unlist(strsplit(basename(tumorBase), "\\."))[1]
copy <- dat$copy

coord <- dat[,c("chr", "start", "end")]
coord$chr <- gsub("chr", "", coord$chr)
coord$chr <- factor(coord$chr, levels=unique(coord$chr))
if(coord$start[1] != 0) {
  coord$start <- coord$start - 1 # start at 0
}

# need to remove NAs, or else the coordinates get all messed up (the start coords sometimes end up greater than end coords if cross chr boundaries)
keepIdx <- !is.na(copy)
copy <- copy[keepIdx]
fullCoord <- coord
coord <- coord[keepIdx,]

CNA.object <- CNA(copy, coord$chr, coord$start, data.type="logratio", sampleid=sampName)

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object)

# remove any duplicates caused by missing data
segment.smoothed.CNA.object$output_withDups <- segment.smoothed.CNA.object$output
segment.smoothed.CNA.object$output <- unique(segment.smoothed.CNA.object$output)
rownames(segment.smoothed.CNA.object$output) <- NULL

# need to convert to bed format, where end coords are exclusive. just add 250000 (assumes always using 250kb windows)
segment.smoothed.CNA.object$output$newEnd <- segment.smoothed.CNA.object$output$loc.end + 250000
segment.smoothed.CNA.object$output$chr <- paste0("chr", segment.smoothed.CNA.object$output$chrom)
write.table(segment.smoothed.CNA.object$output[,c("chr", "loc.start", "newEnd", "seg.mean")], file=paste0(outputBase, ".dnacopy.bed"), quote=F, sep="\t", col.names=F, row.names=F)

