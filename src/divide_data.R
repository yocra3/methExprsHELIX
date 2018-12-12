###############################################################################
#' Divide GenomicRatioSet (methylation) and SummarizedExperiment (Expression) to run
#' linear models
#' @param methy Path with the GenomicRatioSet
#' @param exprs Path with the SummarizedExperiment
#' @param model Model used in the linear model (see Define linear models section)
#' @param out_fold Path with the folder to output the results
#' @example 
#' Rscript divide_datasets.R ./data/gset_sm_res1.RData ./data/summ_exp_sm_res1.RData ./results/model1
###############################################################################

## Load libraries ####
library('minfi')

arg <- commandArgs(trailingOnly = T)

## Load methylation ####
methy <- arg[1]
load(methy) # gset_res

## Load Expression ####
exprs <- arg[2]
load(exprs) # se_res

# Check data consistency ####
stopifnot(colnames(gset_res) == colnames(se_res), "Sample ids in methylation dataset must be equal to sample ids in expression dataset")

out_fold <- arg[3]

print("Data loaded")

## Get data matrices
ma <- assays(gset_res)$Beta # methylation assay
ea <- assays(se_res)$exprs # expression assay

## Get probes ranges
mr <- rowRanges(gset_res)[,0] # methylation ranges
er <- rowRanges(se_res)[,0] # expression ranges

# Automatic subsetting + saving
count = 0
for (x in 1:22) {
  chr <- paste('chr', x, sep = '')
  mrsub <- mr[seqnames(mr) == chr]
  ersub <- er[seqnames(er) == chr]
  fo <- findOverlaps(mrsub + 5e5, ersub)
  
  cpgnames <- names(mrsub)
  tcnames <- names(ersub)
  cpg <- cpgnames[from(fo)]
  tc <- tcnames[to(fo)]
  
  overlaps <- cbind(cpg, tc)
  masub <- ma[cpgnames,]
  easub <- ea[tcnames,]
  save(overlaps, file = paste0(out_fold, "/", 'overlaps', x, '.RData'))
  save(masub, file = paste0(out_fold, "/", 'masub', x, '.RData'))
  save(easub, file = paste0(out_fold, "/", 'easub', x, '.RData'))
  print(paste(x, nrow(overlaps), sep = ': '))
  count <- count + nrow(overlaps)
}
print(count)

########################### Test with CHROMOSOME 22
# mrsub <- mr[seqnames(mr) == 'chr22'] # mr subset
# ersub <- er[seqnames(er) == 'chr22'] # er subset
# fo <- findOverlaps(mrsub + 5e5, ersub) # overlaps
# 
# cpgnames <- names(mrsub) # cpg names of the subset
# tcnames <- names(ersub) # tc names of the subset
# cpg <- cpgnames[from(fo)] # overlapped cpgs
# tc <- tcnames[to(fo)] # overlapped tcs
# 
# overlaps <- cbind(cpg, tc) # 308074
# masub <- ma[cpgnames,] # [1] 6831  832
# easub <- ea[tcnames,] # [1] 1382  832
###################################################
