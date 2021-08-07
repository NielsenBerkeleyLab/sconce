#!/usr/bin/env Rscript
library(HMMcopy)

# Wed 10 Mar 2021 04:56:04 PM PST
# script to run HMMcopy, but without GC and mappability corrections.
# tumor read counts are normalized by matched normal sample (averaged),
# then HMMcopy is run

# results are written to the same directory as the input

args <- commandArgs(trailingOnly=T)

depthFile <- args[1]
healthyFile <- args[2]

tumor <- read.table(depthFile)
colnames(tumor) <- c("chr", "start", "end", "reads")
tumor$copy <- log2(tumor$reads / median(tumor$reads))
tumor$copy[is.infinite(tumor$copy)] <- NA

healthy <- read.table(healthyFile)
colnames(healthy) <- c("chr", "start", "end", "reads")
healthy$copy <- log2(healthy$reads / median(healthy$reads))
healthy$copy[is.infinite(healthy$copy)] <- NA

somatic <- tumor
somatic$copy <- tumor$copy - healthy$copy
somatic_seg <- HMMsegment(somatic)

write.table(somatic, quote=F, col.names=T, row.names=F, sep=",", file=gsub("depth", "csv", depthFile))
somatic <- cbind(somatic, state=somatic_seg$state-1) # subtract 1 because state 1 means <=0 copies, stage 2 means 1 copy, etc
write.table(somatic[,c("chr", "start", "end", "state")], col.names=F, row.names=F, quote=F, sep="\t", file=gsub("depth", "hmmcopyBinnedSegs.bed", depthFile))

