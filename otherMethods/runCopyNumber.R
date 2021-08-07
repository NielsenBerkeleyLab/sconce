library(copynumber)
options(scipen=500) # turn off scientific notation; default is 0

# Tue 20 Apr 2021 04:57:43 PM PDT
# script to collect csv files from HMMcopy and then run copynumber

# Arg 1: path to list of tumor file basenames (auto appends .csv to basenames)
# Arg 2: /path/to/output/csv

args <- commandArgs(trailingOnly=T)
tumorBase <- read.table(args[1], stringsAsFactors=F)$V1
outputBase <- args[2]

coord <- read.table(paste0(tumorBase[1], ".csv"), sep=",", stringsAsFactors=F, header=T)[,c("chr", "start")]
coord$chr <- gsub("chr", "", coord$chr)
if(coord$start[1] != 0) {
  coord$start <- coord$start - 1 # start at 0
}

allTumors <- do.call(cbind, lapply(tumorBase, FUN=function(f) {
  dat <- read.table(paste0(f, ".csv"), sep=",", stringsAsFactors=F, header=T)
  sampName <- unlist(strsplit(basename(f), "\\."))[1]
  copy <- data.frame(dat$copy)
  colnames(copy) <- sampName
  copy
}))

allTumors <- cbind(coord, allTumors)

# impute missing data
imp.data <- imputeMissing(data=allTumors,method="constant")

# winsorize data to handle the outliers
wins.imp <- winsorize(imp.data)

windowSize <- wins.imp$pos[2]
multipcf.segments <- multipcf(data=wins.imp, Y=allTumors, return.est=T, digits=6)

multipcf.segments$estimates$end <- multipcf.segments$estimates$pos + windowSize
colnames(multipcf.segments$estimates)[colnames(multipcf.segments$estimates) == "pos"] <- "start"
multipcf.segments$estimates <- multipcf.segments$estimates[,c(1:2,ncol(multipcf.segments$estimates), 3:(ncol(multipcf.segments$estimates)-1))]

write.csv(cbind(multipcf.segments$estimates[,1:3], 2^multipcf.segments$estimates[,4:ncol(multipcf.segments$estimates)]), file=paste(outputBase,".multi.exp.csv",sep=""), row.names = FALSE, quote=F)

pcf.segments <- cbind(wins.imp[,1:2], do.call(cbind, lapply(colnames(wins.imp)[3:ncol(allTumors)], FUN=function(x) {
  curr <- pcf(data=wins.imp[,c("chrom", "pos", x)], Y=allTumors[,c("chr", "start", x)], return.est=T, normalize=T, digits=6)
  curr$estimates[,x]
})))
colnames(pcf.segments) <- colnames(wins.imp)

pcf.segments$end <- pcf.segments$pos + windowSize
colnames(pcf.segments)[colnames(pcf.segments) == "pos"] <- "start"
pcf.segments <- pcf.segments[,c(1:2,ncol(pcf.segments), 3:(ncol(pcf.segments)-1))]

write.csv(cbind(pcf.segments[,1:3], 2^pcf.segments[,4:ncol(pcf.segments)]), file=paste(outputBase,".indv.exp.csv",sep=""), row.names = FALSE, quote=F)


