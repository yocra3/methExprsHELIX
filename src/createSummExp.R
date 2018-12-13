###############################################################################
#' Create SummarizedExperiment for analysis (Gene Expression)
###############################################################################

library(Biobase)
library(SummarizedExperiment)
library(sva)
library(isva)
library(limma)    
library(SmartSVA) 

## Load SummarizedExperiment
load("./data/transcriptome_subcohort_notfitr_inclsex_v3.RData")
gexp <- transcriptome_subcohort_notfitr_inclsex

## Load common IDs
load("results/preprocessFiles/comIds.Rdata")

# Select common IDs
gexp <- gexp[ , comIds]

## Remove probes with low call rate
gexp <- gexp[fData(gexp)$fil1 == "no", ]

# Making SE
se <- makeSummarizedExperimentFromExpressionSet(gexp)
save(se, file = "results/preprocessFiles/Expression_SE_raw.RData")

## Compute Residuals ####
mat <- assay(se)
pd <- colData(se)

# Protegim les variables que llavors voldrem eliminar a traves dels residuals
Y.r <- t(resid(lm(t(mat) ~ cohort + e3_sex + age_sample_years +
                    NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                    Mono_6, data = pd)))

# Calculem quantes SVs necessitem
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
print(n.sv)
# 57 SVs

# Busquem les SVs
mod <- model.matrix(~ cohort + e3_sex + age_sample_years +
                      NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                      Mono_6, data = pd)
sv.obj <- smartsva.cpp(mat, mod, mod0 = NULL, n.sv = n.sv)

# Eliminem l'efecte de les variables que hem seleccionat i SVs,
# es a dir, obtenim els residuals
model <- cbind(mod, sv.obj$sv)
res <- residuals(lmFit(mat, model), mat)

# Save
assay(se) <- res
save(se, file = "results/preprocessFiles/Expression_SE_residuals.RData")
