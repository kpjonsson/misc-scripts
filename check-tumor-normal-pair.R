#!/usr/bin/env Rscript
# Quick QC of distribution of allele fractions of SNPs in  tumor-normal pair
suppressWarnings(library(data.table))
suppressWarnings(library(Cairo))
args = commandArgs(TRUE)
if (length(args) < 1) stop('Usage: check-TN-pair.R countsMerged_sample.dat.gz')

pileup = read.csv(args[1],
                  stringsAsFactors = FALSE,
                  colClasses = rep(c("character", "numeric","character", "numeric"), c(1,1,2,8)))
name = gsub('^_', '', gsub('.dat.gz', '', gsub('countsMerged_', '', basename(args[1]))))

ii = which(pileup$File1E <= Inf & pileup$File1D <= Inf & pileup$File2E <= Inf & pileup$File2D <= Inf)
rcmat = pileup[ii, 1:2]
rcmat$NOR.DP = pileup$File1R[ii] + pileup$File1A[ii]
rcmat$NOR.RD = pileup$File1R[ii]
rcmat$TUM.DP = pileup$File2R[ii] + pileup$File2A[ii]
rcmat$TUM.RD = pileup$File2R[ii]
rcmat0 = rcmat[rcmat$NOR.DP >= 25,]
rcmat0$vafN = rcmat0$NOR.RD/rcmat0$NOR.DP
rcmat0$vafT = rcmat0$TUM.RD/rcmat0$TUM.DP

ii = which(abs(rcmat0$vafN - 0.5) < 0.45)
CairoPNG(paste0(name, '-TN-check.png'), w = 500, h = 500)
plot(rcmat0$vafN[ii], rcmat0$vafT[ii], pch=".", cex=2)
dev.off()

