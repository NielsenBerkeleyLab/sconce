library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19)

# Mon 26 Apr 2021 04:24:44 PM PDT
# script to run AneuFinder on bam files. should be followed by reformatAneufinderOutput.R to get bed files
# Arg 1: input directory (links to all bam and bai files)
# Arg 2: output dir
# Arg 3: optional flag: should we use chr? default: yes

args <- commandArgs(trailingOnly=T)

datafolder <- args[1]
outputfolder <- args[2]

if(!dir.exists(outputfolder)) {
  dir.create(outputfolder)
}

chromosomes <- paste0("chr", c(1:22,'X','Y'))
# if pass anything as second arg, assume the flag is to not use chr in chromosome names (ie default is to use chr)
if(!is.na(args[3])) {
  chromosomes <- c(1:22, 'X', 'Y')
}
system.time(Aneufinder(inputfolder=datafolder, outputfolder=outputfolder, assembly='hg19', numCPU=14, binsizes=c(250000), chromosomes=chromosomes, correction.method='GC', GC.BSgenome=BSgenome.Hsapiens.UCSC.hg19))
# output is "$f"_aneuOut/BROWSERFILES/method-edivisive/binsize_2e+05_stepsize_2e+05_CNV.bed.gz

