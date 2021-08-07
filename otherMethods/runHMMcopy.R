library(HMMcopy)
options(scipen=500) # turn off scientific notation; default is 0

# Tue 20 Apr 2021 04:14:15 PM PDT
# script to run HMMcopy, following https://www.bioconductor.org/packages/release/bioc/vignettes/HMMcopy/inst/doc/HMMcopy.pdf
# usage:
#   Rscript --vanilla runHMMcopy.R "$wig" "$avgDiploidWig" "$gc" "$map"

args <- commandArgs(trailingOnly=T)

tumorWig <- args[1]
avgDiploidWig <- args[2]
gcFile <- args[3]
mapFile <- args[4]

# read in files
tumorReads <- wigsToRangedData(tumorWig, gcFile, mapFile)
avgDiploidReads <- wigsToRangedData(avgDiploidWig, gcFile, mapFile)

# do read correction
tumorCopy <- correctReadcount(tumorReads)
avgDiploidCopy <- correctReadcount(avgDiploidReads)

# do matchted tumor-normal correction
somaticCopy <- tumorCopy
somaticCopy$copy <- tumorCopy$copy - avgDiploidCopy$copy

# do segmentation
somaticSegments <- HMMsegment(somaticCopy)

# save output
# 0 indexed
if(somaticCopy$start[1] != 0) {
  somaticCopy$start <- somaticCopy$start - 1
  somaticCopy$end <- somaticCopy$end - 1
}
write.table(somaticCopy, quote=F, col.names=T, row.names=F, sep=",", file=gsub("wig", "csv", tumorWig))
somaticCopy <- cbind(somaticCopy, state=somaticSegments$state-1) # subtract 1 because state 1 means <=0 copies, stage 2 means 1 copy, etc
write.table(somaticCopy[,c("chr", "start", "end", "state")], col.names=F, row.names=F, quote=F, sep="\t", file=gsub("wig", "hmmcopyBinnedSegs.bed", tumorWig))

