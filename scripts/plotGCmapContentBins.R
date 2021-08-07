library(ggplot2)
library(reshape2)
library(cowplot)
library(scales)

# Tue 13 Apr 2021 12:34:05 PM PDT
# check if adding gc [and mappability] content makes our predictions of per bin tumor counts any better
# without any library downsampling, for hg19

# compare (via sum of squares)
#   tumor_ij ~ diploidMean_j        (proxy for our method w/o gc content)
#   tumor_ij ~ diploidMean_j + gc_j (add gc content)
#   tumor_ij ~ gc_j                 (just gc content)

# if passed a mappability file, also does:
#   tumor_ij ~ diploidMean_j + mappability      (add mappability content)
#   tumor_ij ~ diploidMean_j + gc + mappability (add gc + mappability content)
#   tumor_ij ~ gc + mappability                 (just gc + mappability content)
#   tumor_ij ~ mappability                      (just mappability content)

# Arg 1: filelist of diploid files
# Arg 2: filelist of tumor files
# Arg 3: /path/to/gc/file, from bedtools nuc (expects gc to be read in to 5th col, named X5_pct_gc)
# Arg 4: /path/to/output/plots
# Arg 5: /path/to/mappability/file, optional. expects a col named "mappability"

# sample usage:
#   Rscript fitGCcontentBins.R diploid_cov_unif_250kb tumor_cov_unif_250kb hg19/hg19_lite.fa.250kb.gc.bed fitGCcontentBins_plots/diploid_cov_unif_250kb_gc hg19/mappability/wgEncodeDukeMapabilityUniqueness35bp.250kb.map.bed > fitGCcontentBins_plots/diploid_cov_unif_250kb_gc.log


args <- commandArgs(trailingOnly=TRUE)
diploidFileList <- read.table(args[1], stringsAsFactors=F)$V1
tumorFileList <- read.table(args[2], stringsAsFactors=F)$V1
gcFile <- args[3]
outputFileBase <- args[4]

# read in all diploid files. Assumes all coords are consistent across all samples
diploid <- do.call(cbind, lapply(diploidFileList, FUN=function(f) {
  tab <- read.table(f, stringsAsFactors=F)
  samp <- unlist(strsplit(basename(f), "\\."))[1]
  colnames(tab) <- c("chr", "start", "end", samp)
  tab[,samp, drop=F]
}))
diploidMean <- rowMeans(diploid)

tumor <- do.call(cbind, lapply(tumorFileList, FUN=function(f) {
  tab <- read.table(f, stringsAsFactors=F)
  samp <- unlist(strsplit(basename(f), "\\."))[1]
  colnames(tab) <- c("chr", "start", "end", samp)
  tab[,samp, drop=F]
}))
tumorNames <- colnames(tumor)

gc <- read.table(gcFile, stringsAsFactors=F, sep="\t", comment.char="", header=T)

# if passed a mappability file
if(length(args) > 4) {
  mappabilityFile <- args[5]
  mappability <- read.table(mappabilityFile, stringsAsFactors=F, sep="\t", header=T)
  diploidGCtumors <- cbind(data.frame(diploidMean=diploidMean, gc=gc[,grepl("_pct_gc", colnames(gc))], mappability=mappability[,"mappability"]), tumor)

  gcCompSS_mat <- sapply(tumorNames, FUN=function(samp) {
    currTab <- diploidGCtumors[,c("gc", "mappability", "diploidMean", samp)]
    diploidOnly <- lm(as.formula(paste(samp, " ~ diploidMean", sep="")), data=currTab)
    diploidGC <- lm(as.formula(paste(samp, " ~ diploidMean + gc", sep="")), data=currTab)
    GCOnly <- lm(as.formula(paste(samp, " ~ gc", sep="")), data=currTab)
    diploidMap <- lm(as.formula(paste(samp, " ~ diploidMean + mappability", sep="")), data=currTab)
    diploidGCmap <- lm(as.formula(paste(samp, " ~ diploidMean + gc + mappability", sep="")), data=currTab)
    gcMap <- lm(as.formula(paste(samp, " ~ gc + mappability", sep="")), data=currTab)
    mapOnly <- lm(as.formula(paste(samp, " ~ mappability", sep="")), data=currTab)
    pred <- data.frame(diploidOnly=predict(diploidOnly, currTab),
                       diploidGC=predict(diploidGC, currTab),
                       GCOnly=predict(GCOnly, currTab),
                       diploidMap=predict(diploidMap, currTab),
                       diploidGCmap=predict(diploidGCmap, currTab),
                       gcMap=predict(gcMap, currTab),
                       mapOnly=predict(mapOnly, currTab))
    diploidOnlySS <- sum((currTab[,samp] - pred$diploidOnly)^2)
    diploidGCSS <- sum((currTab[,samp] - pred$diploidGC)^2)
    gcOnlySS <- sum((currTab[,samp] - pred$GCOnly)^2)
    diploidMapSS <- sum((currTab[,samp] - pred$diploidMap)^2)
    diploidGCmapSS <- sum((currTab[,samp] - pred$diploidGCmap)^2)
    gcMapSS <- sum((currTab[,samp] - pred$gcMap)^2)
    mapOnlySS <- sum((currTab[,samp] - pred$mapOnly)^2)
    c(diploidOnly=diploidOnlySS, diploidGC=diploidGCSS, gcOnly=gcOnlySS, diploidMap=diploidMapSS, diploidGCmap=diploidGCmapSS, gcMap=gcMapSS, mapOnly=mapOnlySS)
  })
} else {
  diploidGCtumors <- cbind(data.frame(diploidMean=diploidMean, gc=gc[,grepl("_pct_gc", colnames(gc))],), tumor)
  gcCompSS_mat <- sapply(tumorNames, FUN=function(samp) {
    currTab <- diploidGCtumors[,c("gc", "diploidMean", samp)]
    diploidOnly <- lm(as.formula(paste(samp, " ~ diploidMean", sep="")), data=currTab)
    diploidGC <- lm(as.formula(paste(samp, " ~ diploidMean + gc", sep="")), data=currTab)
    GCOnly <- lm(as.formula(paste(samp, " ~ gc", sep="")), data=currTab)
    pred <- data.frame(diploidOnly=predict(diploidOnly, currTab),
                       diploidGC=predict(diploidGC, currTab),
                       GCOnly=predict(GCOnly, currTab))
    diploidOnlySS <- sum((currTab[,samp] - pred$diploidOnly)^2)
    diploidGCSS <- sum((currTab[,samp] - pred$diploidGC)^2)
    gcOnlySS <- sum((currTab[,samp] - pred$GCOnly)^2)
    c(diploidOnly=diploidOnlySS, diploidGC=diploidGCSS, gcOnly=gcOnlySS)
  })

}

gcCompSS <- as.data.frame(t(gcCompSS_mat))
gcCompSS$tumor <- rownames(gcCompSS)
toPlot <- melt(gcCompSS)

print(ks.test(gcCompSS$diploidOnly, gcCompSS$diploidGC))
print(ks.test(gcCompSS$diploidOnly, gcCompSS$gcOnly))
print(ks.test(gcCompSS$diploidGC, gcCompSS$gcOnly))

if(length(args) > 4) {
  print(ks.test(gcCompSS$diploidOnly, gcCompSS$diploidMap))
  print(ks.test(gcCompSS$diploidOnly, gcCompSS$diploidGCmap))
  print(ks.test(gcCompSS$diploidOnly, gcCompSS$gcMap))
  print(ks.test(gcCompSS$diploidOnly, gcCompSS$mapOnly))
}


# rename variables for paper plot
colnames(toPlot) <- c("tumor", "model", "value")
toPlot$model <- factor(as.character(toPlot$model), levels=c("diploidOnly", "diploidGC", "diploidMap", "diploidGCmap", "gcOnly", "mapOnly", "gcMap"))

shortLabs <- list(diploidOnly=bquote(A),
             diploidGC=bquote(B),
             diploidMap=bquote(C),
             diploidGCmap=bquote(D),
             gcOnly=bquote(E),
             mapOnly=bquote(F),
             gcMap=bquote(G))

labs <- list(diploidOnly=bquote(A: x[iA] %~% mu[i]),
             diploidGC=bquote(B: x[iA] %~% mu[i] + zeta[i]),
             diploidMap=bquote(C: x[iA] %~% mu[i] + eta[i]),
             diploidGCmap=bquote(D: x[iA] %~% mu[i] + zeta[i] + eta[i]),
             gcOnly=bquote(E: x[iA] %~% zeta[i]),
             mapOnly=bquote(F: x[iA] %~% eta[i]),
             gcMap=bquote(G: x[iA] %~% zeta[i] + eta[i]))


# make into same plot
boxPlot <- ggplot(toPlot, aes(x=model, y=value, colour=model)) + geom_boxplot(notch=T) + guides(colour=F) + labs(x="model", y="sum of squared differences") + theme_bw() + scale_colour_manual(breaks=names(labs), labels=labs, values=hue_pal()(length(labs))) + scale_x_discrete(labels=shortLabs)

ecdfPlot <- ggplot(toPlot, aes(x=value, group=model, colour=model)) + stat_ecdf() + labs(x="sum of squared differences", y="ECDF") + theme_bw() + scale_colour_manual(breaks=names(labs), labels=labs, values=hue_pal()(length(labs)))

legend <- get_legend(ecdfPlot + theme(legend.box.margin=margin(0, 0, 0, 0)))

pGrid <- plot_grid(plotlist=list(boxPlot, ecdfPlot + theme(legend.position="none")), align='vh', labels="AUTO", ncol=2, rel_widths=c(1, 1))
toSave <- plot_grid(pGrid, legend, rel_widths = c(1, .3))
save_plot(paste0(outputFileBase, "sumSqErr_ecdf.png"), toSave)

