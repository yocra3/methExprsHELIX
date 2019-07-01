###############################################################################
# Test cell type specific associations 
###############################################################################
library(dplyr)
library(ggplot2)

load("results/preprocessFiles/gexpAnnotation.Rdata")

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df
featsU <- featStatsDF

load("results/MethComBatExpResidualsCellAdj/allres_simP_cpgs.Rdata")
modC <- df
featsC <- featStatsDF


mergeTB <- modU %>%
  left_join(modC, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  filter(sigPair.x == TRUE | sigPair.y == TRUE) %>%
  mutate(sigType = ifelse(sigPair.x == TRUE, ifelse(sigPair.y == TRUE, "Both", "Crude"), "Cell"))

## Get candidate pairs
arrange(mergeTB, desc(-log10(p.value.x) + log10(p.value.y))) %>% select(CpG, TC, starts_with("FC"), starts_with("p.value"))

## Test: cg02947214 TC14002214.hg.1

## Load data
load("results/MethComBatExpResidualsNoCellAdj/easubchr14.RData")
load("results/MethComBatExpResidualsNoCellAdj/masubchr14.RData")
load("results/MethComBatExpResidualsNoCellAdj/pheno.RData")


plot(masub["cg02947214", ], easub["TC14002214.hg.1", ])


models <- c(cell = easub[tc,] ~ masub[cpg,] + pheno$cohort + pheno$e3_sex +
              pheno$age_sample_years  + pheno$NK_6 + pheno$Bcell_6 +
              pheno$CD4T_6 + pheno$CD8T_6 + pheno$Gran_6 + pheno$Mono_6, 
            nocell = easub[tc,] ~ masub[cpg,] + pheno$cohort + pheno$e3_sex +
              pheno$age_sample_years,
            cellStrat = easub[tc,] ~ masub[cpg,] + pheno$cohort + 
              pheno$age_sample_years  + pheno$NK_6 + pheno$Bcell_6 +
              pheno$CD4T_6 + pheno$CD8T_6 + pheno$Gran_6 + pheno$Mono_6, 
            nocellStrat = easub[tc,] ~ masub[cpg,] + pheno$cohort + 
              pheno$age_sample_years)

cpg <- "cg02947214"
tc <- "TC14002214.hg.1"
summary(lm(models$nocell, data = environment()))
summary(lm(models$cell, data = environment()))


intmod <- easub[tc,] ~ masub[cpg,] + pheno$cohort + pheno$e3_sex +
  pheno$age_sample_years +  masub[cpg,]*pheno$NK_6 + masub[cpg,]*pheno$Bcell_6 +
  masub[cpg,]*pheno$CD4T_6 + masub[cpg,]*pheno$CD8T_6 + masub[cpg,]*pheno$Gran_6 + 
  masub[cpg,]*pheno$Mono_6 - 1
summary(lm(intmod, data = environment()))

df_dat <- data.frame(pheno, tc = easub[tc,], cpg = masub[cpg,])
ggplot(df_dat, aes(x = cpg, y = tc, col = Bcell_6)) + geom_point()
ggplot(df_dat, aes(x = cpg, y = tc, col = CD8T_6 > 0.13)) + geom_point()


cpg <- "cg08766149"
tc <- "TC14000982.hg.1"

cpg <- "cg02947214"
tc <- "TC14000841.hg.1"


## Get candidate pairs (Cell >> Crude)
arrange(mergeTB, desc(-log10(p.value.y) + log10(p.value.x))) %>%
  select(CpG, TC, starts_with("FC"), starts_with("p.value")) %>% 
  filter(grepl("TC14", TC))


cpg <- "cg23090046"
tc <- "TC14002312.hg.1"


intmod <- easub[tc,] ~ masub[cpg,] + pheno$cohort + pheno$e3_sex +
  pheno$age_sample_years +  masub[cpg,]*pheno$CD4T_6 +
  masub[cpg,]*pheno$CD8T_6 + masub[cpg,]*pheno$Gran_6 + 
  - 1
models <- c(cell = easub[tc,] ~ masub[cpg,] + pheno$cohort + pheno$e3_sex +
              pheno$age_sample_years  + pheno$CD4T_6 + pheno$CD8T_6 + 
              pheno$Gran_6, 
            nocell = easub[tc,] ~ masub[cpg,] + pheno$cohort + pheno$e3_sex +
              pheno$age_sample_years)
summary(lm(models$nocell, data = environment()))
summary(lm(models$cell, data = environment()))
summary(lm(intmod, data = environment()))

