###############################################################################
#' Create table with population descriptives
###############################################################################

library(dplyr)
library(S4Vectors)

## Load data  ####
load("results/MethComBatExpResidualsCellAdj/pheno.RData")

getSum <- function(vec, type = "continuous"){
  if (type == "continuous"){
    c(mean = mean(vec, na.rm = TRUE), range = range(vec, na.rm = TRUE))
  } else if (type == "categorical"){
    t <- table(vec)
    data.frame(names = names(t),  tab = as.vector(t), 
               props = as.vector(prop.table(t)))
  }
}

## Categorical variables ####
### Describe Cohort ("cohort") and Sex ("e3_sex")
catVars <- c("cohort", "e3_sex")
catVecs <- lapply(catVars, function(x) pheno[, x])
catSums <- lapply(catVecs, getSum, type = "categorical")

## Create table
catTab <- Reduce(rbind, catSums)
## Convert proportions to percentages
catTab$props <- paste0(sprintf("%.2f", catTab$props*100), "%")
write.table(catTab, file = "paper/PopDescrip_cat.txt", col.names = TRUE, 
            quote = FALSE, row.names = FALSE)

## Categorical variables ####
### Describe Age ("age_sample_years"), zBMI ("hs_zbmi_theano") and 
### Cell types ("NK", "Bcell", "CD4T", "CD8T", "Eos", "Mono", "Neu")
contVars <- c("age_sample_years", "hs_zbmi_theano", "NK", "Bcell", "CD4T", "CD8T", "Eos", "Mono", "Neu")
contVecs <- lapply(contVars, function(x) pheno[, x])
contSums <- lapply(contVecs, getSum, type = "continuous")

## Create table
contTab <- data.frame(Reduce(rbind, contSums))
## Convert proportions to percentages
contVals <- sprintf("%.3f (%.3f-%.3f)", contTab$mean, contTab$range1, contTab$range2)
contTab <- data.frame(names = contVars, vals = contVals)
write.table(contTab, file = "paper/PopDescrip_cont.txt", col.names = TRUE, 
            quote = FALSE, row.names = FALSE)
