library(ggplot2)
library(reshape2)

# Thu 12 Aug 2021 12:20:06 PM PDT
# script to plot one genome trace from SCONCE output
#
# Arg 1: /path/to/healthy/average/bed/file
# Arg 2: /path/to/observed/tumor/read/depths/bed/file
# Arg 3: /path/to/sconce/output/bed/file
# Arg 4: /path/to/output/plot
# Arg 5: quoted text for the plot title (ex sample name)
# Arg 6: [optional] /path/to/ground/truth/ploidy/bed/file (ie for simulations)
#
# sample usage:
#   Rscript scripts/plotGenomeTrace.R test/test_healthy_avg.bed test/test_cancer_cell.bed test/ref_output_k5.bed test/ref_plot_k5.png "Reference Genome Trace for SCONCE (k=5)" test/true_cancer_cell.bed

# read in command line args
args <- commandArgs(trailingOnly=TRUE)
avgDiploidFile <- args[1]
obsReadDepthFile <- args[2]
sconceFile <- args[3]
plotFile <- args[4]
plotTitle <- args[5]

# check if given a simulated ground truth file
groundTruthFile <- NA
if(!is.na(args[6]) && nchar(args[6]) != 0) {
  groundTruthFile <- args[6]
}

# read in files
bedHeader <- c("chr", "start", "end", "estPloidy")
obsReadDepth <- read.table(obsReadDepthFile, stringsAsFactors=F)
colnames(obsReadDepth) <- c("chr", "start", "end", "readDepth")
obsReadDepth$idx <- 1:nrow(obsReadDepth)
coordOrdering <- obsReadDepth[,c("chr", "start", "idx")]

avgDiploid <- read.table(avgDiploidFile, stringsAsFactors=F)
colnames(avgDiploid) <- c("chr", "start", "end", "readDepth", "var")
avgDiploid <- merge(coordOrdering, avgDiploid, by=c("chr", "start"))

sconceDat <- read.table(sconceFile, stringsAsFactors=F, header=F)
colnames(sconceDat) <- bedHeader
sconceDat <- merge(coordOrdering, sconceDat, by=c("chr", "start"), all=T)
if(nrow(sconceDat) != nrow(obsReadDepth)) {
  warning(paste0("Could not merge on coordinates for ", sconceFile))
}

sumSq <- NA
if(!is.na(groundTruthFile)) {
  groundTruth <- read.table(groundTruthFile, stringsAsFactors=F)
  colnames(groundTruth) <- c("chr", "start", "end", "avgPloidy")
  groundTruth <- merge(coordOrdering, groundTruth, by=c("chr", "start"))

  # calc sumSqDiff if have groundTruth
  mergedTruthSconce <- merge(groundTruth, sconceDat, by=c("chr", "start"))
  sumSq <- sprintf("SSD = %.2f", round(sum((mergedTruthSconce$avgPloidy - mergedTruthSconce$estPloidy)^2), digits=2))
}

# calc chr boundaries. midpoint calculations done on the fly below, from https://stackoverflow.com/a/54147509
chrBoundaryIdx <- c(1, cumsum(rle(coordOrdering$chr)$lengths))

# scale ploidy into observed read depth range for plotting purposes
scalingFactor <- sum(obsReadDepth$readDepth) / mean(sconceDat[,"estPloidy"], na.rm=T) / nrow(obsReadDepth)
if(!is.na(groundTruthFile)) {
  scalingFactor <- sum(obsReadDepth$readDepth) / mean(groundTruth$avgPloidy) / nrow(obsReadDepth) # for true expected read count
}

# construct plot
trace <- ggplot() +
  geom_line(dat=avgDiploid, aes(x=idx, y=readDepth), colour="lightblue", alpha=0.4) +
  geom_ribbon(dat=avgDiploid, aes(x=idx, ymin=(readDepth - sqrt(var)), ymax=(readDepth + sqrt(var))), fill="lightblue", alpha=0.2) +
  geom_point(dat=obsReadDepth, aes(x=idx, y=readDepth), alpha=0.3, colour="darkgray", size=0.5) +
  geom_line(data=sconceDat, aes(x=idx, y=estPloidy * scalingFactor), colour="#C49A00") +
  scale_y_continuous(sec.axis=sec_axis(~ . / scalingFactor, name="ploidy")) +
  theme_bw() +
  geom_vline(xintercept=chrBoundaryIdx, colour="gray30") +
  theme(legend.position="none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x=element_blank()) +
  scale_x_continuous(breaks=c(chrBoundaryIdx[-length(chrBoundaryIdx)] + diff(chrBoundaryIdx)/2)[c(-19,-21)], labels=c(1:18,20,22, "X", "Y"), expand=c(0,0)) +
  labs(title=plotTitle)
if(!is.na(groundTruthFile)) {
  trace <- trace + geom_line(data=groundTruth, aes(x=idx, y=avgPloidy * scalingFactor), colour="red", linetype="dashed") # add ground truth line
  trace <- trace + labs(subtitle=sumSq, x=NULL, y="read depth")
} else {
  trace <- trace + labs(x=NULL, y="read depth")
}

ggsave(trace, file=plotFile, width=8, height=3)

