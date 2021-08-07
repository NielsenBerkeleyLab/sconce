library(AneuFinder)
library(genomation)

# Mon 26 Apr 2021 08:42:21 PM PDT
# script to run AneuFinder without any GC or mappability correction (ie for simulations)
# Arg 1: input directory (links to all bam and bai files)
# Arg 2: output dir

args <- commandArgs(trailingOnly=T)
aneuInDir <- args[1]
aneuOutDir <- args[2]

fileList <- system(paste0("find ", aneuInDir, " -name \"simu_cancer_cell*hg19_lite.depth\" | sort -V"), intern=T)

if(!dir.exists(aneuOutDir)) {
  dir.create(aneuOutDir)
}

sapply(fileList, FUN=function(f) {
  cellName <- gsub(".depth", "", basename(f))
  bed <- readGeneric(f, chr=1, start=2, end=3, strand=NULL, meta.cols=list(counts=4))
  rDataFile <- paste0(aneuOutDir, "/", cellName, ".aneu.RData")
  save(bed, file=rDataFile)
  model <- findCNVs(rDataFile, method = "edivisive", R=10, sig.lvl=0.1)
  cnvBed <- data.frame(chr=seqnames(model$segments), start=start(model$segments), end=end(model$segments), cnv=elementMetadata(model$segments)$copy.number)
  cnvBedName <- paste0(aneuOutDir, "/", cellName, ".aneu.bed")
  write.table(cnvBed, cnvBedName, sep="\t", quote=F, col.names=F, row.names=F)
})

