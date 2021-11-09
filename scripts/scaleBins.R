# script to scale and shift relative copy number estimates results (ex from CopyNumber and DNAcopy)
# on simulated data by minimizing the sum of squared differences
# so they can be compared to other methods

# Arg 1: unscaled bed file
# Arg 2: ground truth file
# Arg 3: output file

args <- commandArgs(trailingOnly=T)
unscaledFile <- args[1]
groundTruthFile <- args[2]
outputFile <- args[3]

unscaled <- read.table(unscaledFile, stringsAsFactors=F)
groundTruth <- read.table(groundTruthFile, stringsAsFactors=F)
colnames(groundTruth) <- c("chr", "start", "end", "avgPloidy")
colnames(unscaled) <- c("chr", "start", "end", "estPloidy")

merged <- merge(groundTruth, unscaled, by=c("chr", "start"), sort=F)

temp <- merged[!is.na(merged$estPloidy),]
scaleShiftFun <- function(params) {
  scale <- params[1]
  shift <- params[2]
  if(any(temp$estPloidy * scale + shift < 0, na.rm=T)) {
    return(10000000)
  }
  sum((temp$avgPloidy - (temp$estPloidy * scale + shift))^2, na.rm=T)
}
optScaling <- optim(par=c(1,10), method="BFGS", fn=scaleShiftFun)
scaled <- merged
scaled$estPloidy <- merged$estPloidy * optScaling$par[1] + optScaling$par[2]

write.table(scaled[,c("chr", "start", "end.x", "estPloidy")], file=outputFile, quote=F, sep="\t", col.names=F, row.names=F)

