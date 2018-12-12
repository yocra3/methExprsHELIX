#####################################################################
# residuals metilacio
# residuals metilacio = sex, age, cohort, cell types (data=combat)

# residuals expression 
## 1r. calcular variables surrogades (SVs)
## 2n. obtenir residuals: 
# residuals expression = sex, age, cohort, cell types, SV

#####################################################################
### Libraries

library(minfi)
library(sva)
library(isva)
library(limma)    # We use lmFit to fit the lineal model
library(SmartSVA) # We want to compute the SVA to correct methylation data
library(Biobase)

setwd("/scratch/smari")

#####################################################################
#### EXPRESSION

load("summ_exp_sm_1.RData") # se
dd <- se
m <- assays(dd)$exprs
pd <- colData(dd)

# Protegim les variables que llavors voldrem eliminar a traves dels residuals
Y.r <- t(resid(lm(t(m) ~ cohort + e3_sex + age_sample_years +
                    NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                    Mono_6, data = pd)))

# Calculem quantes SVs necessitem
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
# 57 SVs

# Busquem les SVs
mod <- model.matrix(~ cohort + e3_sex + age_sample_years +
                      NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                      Mono_6, data = pd)
sv.obj <- smartsva.cpp(m, mod, mod0 = NULL, n.sv = n.sv)

# Eliminem l'efecte de les variables que hem seleccionat i SVs,
# es a dir, obtenim els residuals
model <- cbind(mod, sv.obj$sv)
res <- residuals(lmFit(m, model), m)

# Save
se_res <- se
assays(se_res)$exprs <- res
save(se_res, file = "summ_exp_sm_res1.RData")

#####################################################################
### METHYLATION

load("gset_sm_1.RData") # gset
dd <- gset
m <- assays(dd)$Beta
pd <- colData(dd)

mod <- model.matrix(~ cohort + e3_sex + age_sample_years +
                      NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                      Mono_6, data = pd)

res <- residuals(lmFit(m, mod), m)

gset_res <- gset
assays(gset_res)$Beta <- res
save(gset_res, file = "gset_sm_res1.RData")


#####################################################################
### EXPRESSION (2)

# load("summ_exp_sm_1.RData")
dd <- se
m <- assays(dd)$exprs
pd <- colData(dd)

# SVs
Y.r <- t(resid(lm(t(m) ~ cohort + e3_sex + age_sample_years + 
                    NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                    Mono_6, data = pd)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
# 57 SVs

mod <- model.matrix(~ cohort + e3_sex + age_sample_years +
                      NK_6 + Bcell_6 + CD4T_6 + CD8T_6 + Gran_6 +
                      Mono_6, data = pd)

sv.obj <- smartsva.cpp(m, mod, mod0 = NULL, n.sv = n.sv)
save(sv.obj, file = "sv.obj_exp_res2.RData")
# load("sv.obj_exp_res2.RData")

res <- residuals(lmFit(m, sv.obj$sv), m)

se_res <- se
assays(se_res)$exprs <- res
save(se_res, file = "summ_exp_sm_res2.RData")

### Association between residuals and covariates

sv.obj$n.sv # 57
sv_s <- sv.obj$sv

pheno <- colData(se)
names(pheno)

phenof <- pheno[,c("cohort", "Period", "e3_sex","age_sample_years",
                   "h_native", # "h_ethnicity_c", "h_ethnicity_cauc",
                   "ethn_dummy","ethn_dummy2","ethn_PC1", "ethn_PC2",
                   # "round_gexp",
                   "r_b","cohort_GE","dil","evap",
                   "date_start",
                   #"round_rna_extr",
                   "extr_batch",
                   "bio_conc_denat_ngul", "bio_rin_denat",
                   "extr_technician", "x260_230","x260_280",
                   "NK", "Bcell", "CD4T", "CD8T", "Eos", "Mono",
                   "Neu", "NK_6", "Bcell_6", "CD4T_6", "CD8T_6",
                   "Mono_6", "Gran_6")]

rsquaresSV <- sapply(1:57, function(r) sapply(phenof, function(y)
  summary(lm(sv_s[, r] ~ y))$adj.r.squared))

write.table(rsquaresSV, file = "rsquaresSV_exp_res2.txt", sep="\t", quote=F)


#####################################################################
### METHYLATION (2)
# Ajustar directament per sexe, cohort, edat i cell type a la
# regressio lineal

