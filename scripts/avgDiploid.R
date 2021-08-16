# Fri 02 Apr 2021 05:09:13 PM PDT
# script to average diploid files together.
# Assumes all files are sorted in the same way (ie doesn't check for coordinate matching)

# Arg 1: file of list of paths to diploid files
# Arg 2: output filename (ex diploid/diploid_avg_cov_unif_200kb.bed)
#
# sample usage:
#   Rscript scripts/fitMeanVarRlnshp.R test/diploidFileList test/test_healthy_avg.bed

args <- commandArgs(trailingOnly=TRUE)

fileList <- read.table(args[1], stringsAsFactors=F)$V1
outputFile <- args[2]

coords <- read.table(fileList[1], stringsAsFactors=F)
colnames(coords) <- c("chr", "start", "end", "depth")
coords[,"depth"] <- NULL

dat <- do.call(cbind, lapply(fileList, FUN=function(f) {
  tab <- read.table(f, stringsAsFactors=F)
  colnames(tab) <- c("chr", "start", "end", f)
  tab[,f, drop=F]
}))

dat2 <- dat
dat2$mean <- apply(dat, 1, mean)
dat2$var <- apply(dat, 1, var)

toWrite <- cbind(coords, dat2$mean, dat2$var)
write.table(toWrite, file=outputFile, col.names=F, row.names=F, quote=F, sep="\t")

