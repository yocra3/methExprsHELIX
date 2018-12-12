###############################################################################
#' Create SummarizedExperiment for analysis (Gene Expression)
###############################################################################

library(Biobase)
library(SummarizedExperiment)

setwd("/Lacie_CRW10023/HELIX_analyses/expr_met_SM")
load("anno_expr_sm_2.RData")
names(anno)
dim(anno)

load("/Lacie_CRW10023/HELIX_preproc/gene_expression/Final_data/transcriptome_subcohort_f1_residuals3_v2.Rdata")
eset <- transcriptome_subcohort_f1_residuals3
feat <- featureData(eset)
eanno <- pData(feat)
names(eanno)
dim(eanno)

# Filter IDs
load("ids.RData")
eset <- eset[,ids] #832

# Filter TCs
# eanno: original expression set annotation
# anno: filtered annotation (chr gl, chr M and haplotypes)
eset <- eset[rownames(anno),] #59130

# Adding anno
pData(feat) <- anno
featureData(eset) <- feat

# Making SE
se <- makeSummarizedExperimentFromExpressionSet(eset)
save(se, file = "summ_exp_sm_1.RData")

