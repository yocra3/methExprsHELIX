###############################################################################
#' Divide GenomicRatioSet (methylation) and SummarizedExperiment (Expression) to run
#' linear models
#' @param methy Path with the GenomicRatioSet
#' @param exprs Path with the SummarizedExperiment
#' @param out_fold Path with the folder to output the results
#' @param Operation Operations to be applied to data (autosomes: select probes in 
#' autosome chromosome; sex-stratify: divide the datasets in males and females)
#' @example 
#' Rscript createInputFiles.R ./data/gset_sm_res1.RData ./data/summ_exp_sm_res1.RData ./results/model1
###############################################################################

## Load libraries ####
library('minfi', quietly = TRUE)
library(SummarizedExperiment, quietly = TRUE)

arg <- commandArgs(trailingOnly = T)

## Load methylation ####
methy <- arg[1]
message("methy File: ", methy)
load(methy) # gset

## Load Expression ####
exprs <- arg[2]
message("Gexp File: ", exprs)
load(exprs) # se

# Check data consistency ####
stopifnot(colnames(gset) == colnames(se), "Sample ids in methylation dataset must be equal to sample ids in expression dataset")

out_fold <- arg[3]

print("Data loaded")

## Remove probes in sexual chromosomes
if ("autosomes" %in% args){
  gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
  se <- se[!seqnames(rowRanges(se)) %in% c("chrX", "chrY"), ]
}

## Remove probes in sexual chromosomes
if ("sex-stratify" %in% args){
  gsetall <- gset
  gset <- gsetall[, gsetall$e3_sex == "male"]
  save(gset, file = paste0(out_fold, "/methyInputMale.Rdata"))

  gset <- gsetall[, gsetall$e3_sex == "female"]
  save(gset, file = paste0(out_fold, "/methyInputFemale.Rdata"))
  
  seall <- se
  se <- seall[, seall$e3_sex == "male"]  
  save(se, file = paste0(out_fold, "/gexpInputMale.Rdata"))

  se <- seall[, seall$e3_sex == "female"]
  save(se, file = paste0(out_fold, "/gexpInputFemale.Rdata"))
} else {
  save(gset, file = paste0(out_fold, "/methyInput.Rdata"))
  save(se, file = paste0(out_fold, "/gexpInput.Rdata"))
}