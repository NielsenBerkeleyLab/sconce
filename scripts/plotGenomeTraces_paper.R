library(ggplot2)
library(reshape2)
library(cowplot)
library(scales)

# Wed 28 Apr 2021 08:27:32 PM PDT
# script to plot genome traces
# Rscript --vanilla "$scripts""/plotGenomeTraces.R" "$datasetName" "$cellID" "$obsReadDepth" "$avgDiploid" "$plotDir" "$sconceFile_k5" "$sconceFile_k8" "$sconceFile_k10" "$hmmCopyFile" "$copyNumberMultiFile" "$copyNumberIndvFile" "$aneuFile" "$groundTruth"

args <- commandArgs(trailingOnly=TRUE)
datasetName <- args[1]
sconceFileKey <- args [2]
cellID <- args[3]
obsReadDepthFile <- args[4]
avgDiploidFile <- args[5]
plotDir <- args[6]
sconceFile_k5 <- args[7]
sconceFile_k10 <- args[8]
sconceFile_k15 <- args[9]
hmmCopyFile <- args[10]
copyNumberMultiFile <- args[11]
copyNumberIndvFile <- args[12]
aneuFile <- args[13]

# if given a ground truth file for sims
if(!is.na(args[14])) {
  groundTruthFile <- args[14]
  if(nchar(args[14]) == 0) {
    groundTruthFile <- NA
  }
}

compFiles <- c(sconceFile_k5, sconceFile_k10, sconceFile_k15, aneuFile, hmmCopyFile, copyNumberMultiFile, copyNumberIndvFile)
compFiles <- compFiles[sapply(compFiles, nchar) > 0]

bedHeader <- c("chr", "start", "end", "estPloidy")

obsReadDepth <- read.table(obsReadDepthFile, stringsAsFactors=F)
colnames(obsReadDepth) <- c("chr", "start", "end", "readDepth")
obsReadDepth$idx <- 1:nrow(obsReadDepth)
coordOrdering <- obsReadDepth[,c("chr", "start", "idx")]

avgDiploid <- read.table(avgDiploidFile, stringsAsFactors=F)
colnames(avgDiploid) <- c("chr", "start", "end", "readDepth", "var")
avgDiploid <- merge(coordOrdering, avgDiploid, by=c("chr", "start"))

if(!is.na(groundTruthFile)) {
  groundTruth <- read.table(groundTruthFile, stringsAsFactors=F)
  colnames(groundTruth) <- c("chr", "start", "end", "avgPloidy")
  groundTruth <- merge(coordOrdering, groundTruth, by=c("chr", "start"))
}

cleanProgramName <- function(f) {
  if(grepl("k5.viterbiDecoded", f)) {
    return("SCONCE (k = 5)")
  } else if(grepl("k10.viterbiDecoded", f)) {
    return("SCONCE (k = 10)")
  } else if(grepl("k15.viterbiDecoded", f)) {
    return("SCONCE (k = 15)")
  } else if(grepl("hmmcopyBinnedSegs", f)) {
    return("HMMcopy")
  } else if(grepl("copynumberBinned.multi", f)) {
    return("CopyNumber (multi)")
  } else if(grepl("copynumberBinned.indv", f)) {
    return("CopyNumber (indv)")
  } else if(grepl("aneufinderBins", f)) {
    return("AneuFinder")
  }
  return(f)
}

dat <- do.call(rbind, lapply(compFiles, FUN=function(f) {
  currEst <- read.table(f, stringsAsFactors=F, header=F)
  colnames(currEst) <- bedHeader

  # merge on start coordinate
  merged <- merge(coordOrdering, currEst, by=c("chr", "start"), all=T)
  if(nrow(merged) != nrow(obsReadDepth)) {
    if(grepl("aneu", f)) {
      print("attempting to subtract 1 from segments that end in 1")
      currEst$start[currEst$start %% 10 == 1] <- currEst$start[currEst$start %% 10 == 1] - 1
      merged <- merge(coordOrdering, currEst, by=c("chr", "start"), all=T)
      if(nrow(merged) != nrow(obsReadDepth)) {
        warning(paste0("Could not merge on coordinates for ", f))
      } else {
        print("successful merge")
      }
    }
  }

  # need to scale copynumber results, since copynumber doesn't output absolute copy number estimates
  if(grepl("copynumber", f, ignore.case=T) && !is.na(groundTruthFile)) {
    temp <- merge(merged, groundTruth)
    optScaling <- optim(par=1, method="Brent", lower=0, upper=1e20, fn=function(scaling) {
      sum((temp$avgPloidy - temp$estPloidy * scaling)^2)
    })
    merged$estPloidy <- merged$estPloidy * optScaling$par
  }

  program <- cleanProgramName(f)
  merged$program <- program

  merged
}))

cnvScaling <- median(obsReadDepth$readDepth) / 2
if(!is.na(groundTruthFile)) {
  cnvScaling <- sum(obsReadDepth$readDepth) / mean(groundTruth$avgPloidy) / nrow(obsReadDepth) # for true expected read count
}

# calc sumSqDiff if have groundTruth
sumSq <- NULL
if(!is.na(groundTruthFile)) {
  sumSq <- do.call(rbind, lapply(unique(dat$program), FUN=function(x) {
    merged <- merge(groundTruth, subset(dat, program == x))
    #data.frame(program=x, sumSq=paste0("sumSq = ", sum((merged$avgPloidy - merged$estPloidy)^2)))
    data.frame(program=x, sumSq=sprintf("(SSD = %.2f)", round(sum((merged$avgPloidy - merged$estPloidy)^2), digits=2)))
  }))
}

# calc chr boundaries
chrBoundaryIdx <- c(1, cumsum(rle(coordOrdering$chr)$lengths))
# midpoint calculations done on the fly below, from https://stackoverflow.com/a/54147509

scalingFactors <-  sapply(unique(dat$program), FUN=function(x) { sum(obsReadDepth$readDepth) / mean(dat[dat$program == x, "estPloidy"], na.rm=T) / nrow(obsReadDepth)})
if(!is.na(groundTruthFile)) {
  cnvScaling <- sum(obsReadDepth$readDepth) / mean(groundTruth$avgPloidy) / nrow(obsReadDepth) # for true expected read count
  scalingFactors <- rep(cnvScaling, length(unique(dat$program)))
  names(scalingFactors) <- unique(dat$program)
}
progColors <- hue_pal()(length(unique(dat$program)))
names(progColors) <- unique(dat$program)

plotList <- lapply(unique(dat$program), FUN=function(currentProgram) {
  currPlot <- ggplot() +
    geom_line(dat=avgDiploid, aes(x=idx, y=readDepth), colour="lightblue", alpha=0.4) +
    geom_ribbon(dat=avgDiploid, aes(x=idx, ymin=(readDepth - sqrt(var)), ymax=(readDepth + sqrt(var))), fill="lightblue", alpha=0.2) +
    geom_point(dat=obsReadDepth, aes(x=idx, y=readDepth), alpha=0.3, colour="darkgray", size=0.5) +
    geom_line(data=dat[dat$program == currentProgram,], aes(x=idx, y=estPloidy * scalingFactors[currentProgram]), colour=progColors[currentProgram]) +
    scale_y_continuous(sec.axis=sec_axis(~ . / scalingFactors[currentProgram], name="ploidy")) +
    theme_bw() +
    geom_vline(xintercept=chrBoundaryIdx, colour="gray30") +
    theme(legend.position="none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x=element_blank()) +
    #scale_x_continuous(breaks=chrBoundaryIdx[-length(chrBoundaryIdx)] + diff(chrBoundaryIdx)/2, labels=c(1:22, "X", "Y"), expand=c(0,0))
    scale_x_continuous(breaks=c(chrBoundaryIdx[-length(chrBoundaryIdx)] + diff(chrBoundaryIdx)/2)[c(-19,-21)], labels=c(1:18,20,22, "X", "Y"), expand=c(0,0))
  if(!is.na(groundTruthFile)) {
    currPlot <- currPlot + geom_line(data=groundTruth, aes(x=idx, y=avgPloidy * scalingFactors[currentProgram]), colour="red", linetype="dashed") # add ground truth line
    #currPlot <- currPlot + geom_text(data=sumSq[sumSq$program == currentProgram,], aes(x=Inf, y=Inf, label=sumSq, vjust=1, hjust=1))
    currPlot <- currPlot + labs(subtitle=paste0(currentProgram, " ", sumSq[sumSq$program == currentProgram,"sumSq"]), x=NULL, y="read depth")
  } else {
    currPlot <- currPlot + labs(subtitle=currentProgram, x=NULL, y="read depth")
  }
  currPlot
})

traces <- plot_grid(plotlist=plotList, ncol=1, align="hv")
if(length(plotList) == 7) {
  save_plot(paste0(plotDir, "/", datasetName, "_", sconceFileKey, "_", cellID, "_genomeTrace.png"), traces, base_width=8, base_height=8)
} else {
  save_plot(paste0(plotDir, "/", datasetName, "_", sconceFileKey, "_", cellID, "_genomeTrace.png"), traces, base_width=8, base_height=8/7*5)
}

