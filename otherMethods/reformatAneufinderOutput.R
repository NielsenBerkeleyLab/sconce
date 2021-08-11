library(AneuFinder)

# Mon 26 Apr 2021 04:25:12 PM PDT
# script to reformat output (from .Rdata files) into bed files
# Arg 1: output directory from runAneufinder.R

args <- commandArgs(trailingOnly=T)
aneuOutDir <- args[1]
fileList <- system(paste0("find ", aneuOutDir, "/MODELS/method-edivisive -name \"*RData\" | sort -V"), intern=T)

cnvDir <- paste0(aneuOutDir, "/cnvCalls")
if(!dir.exists(cnvDir)) {
  dir.create(cnvDir)
}
lapply(fileList, FUN=function(f) {
  load(f)
  currOut <- paste0(cnvDir, "/", gsub(".RData", ".bed", basename(f)))
  bed <- data.frame(chr=seqnames(model$segments), start=start(model$segments), end=end(model$segments), cnv=elementMetadata(model$segments)$copy.number)
  write.table(bed, currOut, sep="\t", quote=F, col.names=F, row.names=F)
})

