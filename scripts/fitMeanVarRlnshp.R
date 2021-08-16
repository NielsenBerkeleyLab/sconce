#!/usr/bin/env Rscript
# Tue 30 Mar 2021 05:00:18 PM PDT
# script to learn the mean and variance relationship coefficients from a set of diploid data, by maximizing the loglikelihood of the negative binomial using the mean/variance 2nd degree polynomial relationship

# Arg 1: file of list of paths to diploid files
# Arg 2: output file 
#
# sample usage:
#   Rscript scripts/fitMeanVarRlnshp.R test/diploidFileList test/test.meanVar

args <- commandArgs(trailingOnly=TRUE)
fileList <- read.table(args[1], stringsAsFactors=F)$V1
outputFile <- args[2]

# read in all diploid files. Assumes all coords are consistent across all samples
dat <- do.call(cbind, lapply(fileList, FUN=function(f) {
  tab <- read.table(f, stringsAsFactors=F)
  colnames(tab) <- c("chr", "start", "end", f)
  tab[,f, drop=F]
}))

totalReads <- sum(dat)
windowProp <- rowSums(dat) / totalReads # window proportions
cellProp <- colSums(dat) / totalReads # cell proportions

expMat <- outer(windowProp, cellProp * totalReads)
rownames(expMat) <- rownames(dat)

loglikeFunMax <- function(param) {
 loglik <- 0
 for(i in 1:length(windowProp)){
   for(j in 1:length(cellProp)) {
     #Eij <- windowProp[i] * cellProp[j] * totalReads
     Eij <- expMat[i, j]
     Vij <- param[1] + param[2] * Eij + param[3] * Eij^2
     p <- (Eij / Vij)
     r <- (Eij^2 /(Vij - Eij))
     if(Eij == 0) {
       if(dat[i,j] == 0) {
         loglik <- loglik + log(1)
       }
       else {
         loglik <- loglik + -1e300
       }
     }
     else {
       if(p < 0 || r < 0) {
         loglik <- loglik + -1e300
       }
       else {
         prob <- dnbinom(x=dat[i,j], size=r, prob=p, log=T)
         if(is.nan(prob)) {
           print(paste0(Eij, ", ", Vij, ", ", p, ", ", r, ", ", dnbinom(x=dat[i,j], size=r, prob=p, log=T), ", [", param[1], ", ", param[2], ", ", param[3], "]"))
           loglik <- loglik + -1e300
         }
         else {
           loglik <- loglik + prob
         }
       }
     }
   }
 }
 return(loglik)
}

print(system.time(res <- optim(c(10,1,1),loglikeFunMax, control=list(fnscale=-1)))) # ~12 mins for Navin 2011 dataset (~60 cells)
sprintf("%.20f", res$par)
meanVarTab <- data.frame(c("intercept", "slope", "poly2"), res$par)
write.table(format(meanVarTab, digits=10), file=outputFile, quote=F, sep="=", col.names=F, row.names=F)

# > res
# $par
# [1] 10.46385712  2.42601322  0.01114518
# 
# $value
# [1] -1847413
# 
# $counts
# function gradient
#      122       NA
# 
# $convergence
# [1] 0
# 
# $message
# NULL

# > sprintf("%.20f", res$par)
# [1] "10.46385711652957084539" "2.42601321762369614987"
# [3] "0.01114518215725581601"

