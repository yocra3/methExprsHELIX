#'##############################################################################
#' Genetic variants and eQTMs
#' This script contains the code for
#'  - Heritability
#'  - Association with genetics
#'##############################################################################

# Load data and libraries ####
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(GenomicRanges)
library(snpStats)
library(MASS)
library(topGO)

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/gexpAnnotation.Rdata")

# Get useful variables ####
CpGsNum <- df %>%
  group_by(CpG) %>%
  summarise(nTC = sum(sigPair)) %>%
  mutate(nCat = ifelse(nTC > 3, "4+", nTC))

codingTCs <- subset(expAnnot, Coding == "coding")$transcript_cluster_id

## Heritability ####
## Create summary tibble
h2df <- read.csv("data/HeritabilityCpGs", as.is = TRUE)
herm <- h2df %>%
  as_tibble() %>%
  ## Remove CpGs with NA in h2
  filter(!is.na(h2_total)) %>%
  mutate(CpG = cgid) %>%
  inner_join(CpGsNum, h2df, by = "CpG") 

h2tot <- herm %>%
  ggplot(aes(x = nCat, y = h2_total)) + 
  geom_boxplot() +
  geom_hline(yintercept = c(0.2, 0.5), linetype="dashed", colour = "blue") + 
  ggtitle("Total Heritability") +
  scale_x_discrete(name = "num of TCs associated with a CpG") +
  scale_y_continuous(name = "h2", limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

h2SNP <- herm %>%
  ggplot(aes(x = nCat, y = h2_total)) + 
  geom_boxplot() + 
  geom_hline(yintercept = c(0.2, 0.5), linetype="dashed", colour = "blue") + 
  ggtitle("SNP Heritability") +
  scale_x_discrete(name = "num of TCs associated with a CpG") +
  scale_y_continuous(name = "h2", limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")



herm %>%
  mutate(sig = ifelse(Combined == "Non-significant", "Non-significant", "Significant"),
         sig = factor(sig, levels = c("Non-significant", "Significant"))) %>%
  glm(h2_total ~ sig, family = "binomial", .) %>%
  summary()
# Call:                                                                                                                                                                   [15/409]
# glm(formula = h2_total ~ sig, family = "binomial", data = .)
# Deviance Residuals:
#   Min       1Q   Median       3Q      Max
# -1.1210  -0.4158  -0.1725   0.2455   1.8392
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)    -1.574966   0.004528 -347.86   <2e-16 ***
#   sigSignificant  1.440844   0.014891   96.76   <2e-16 ***
#   
herm %>%
  mutate(sig = ifelse(Combined == "Non-significant", "Non-significant", "Significant"),
         sig = factor(sig, levels = c("Non-significant", "Significant"))) %>%
  glm(h2_SNPs ~ sig, family = "binomial", .) %>%
  summary()
# Call:
#   glm(formula = h2_SNPs ~ sig, family = "binomial", data = .)
# 
# Deviance Residuals:
#   Min        1Q    Median        3Q       Max
# -0.70686  -0.34995  -0.30398   0.09453   2.31096
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)    -2.762317   0.007221 -382.56   <2e-16 ***
#   sigSignificant  1.502830   0.018521   81.14   <2e-16 ***
#   


herm %>%
  filter(Type != "Non-significant") %>%
  glm(h2_total ~ Type, family = "binomial", .) %>%
  summary()
# Call:
#   glm(formula = h2_total ~ Type, family = "binomial", data = .)
# Deviance Residuals:
#   Min        1Q    Median        3Q       Max
# -1.23705  -0.34753  -0.01448   0.34700   1.22204
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.29345    0.01797  -16.33   <2e-16 ***
#   TypeMulti    0.43261    0.02955   14.64   <2e-16 ***
#   
herm %>%
  filter(Type != "Non-significant") %>%
  glm(h2_SNPs ~ Type, family = "binomial", .) %>%
  summary()
# Call:
#   glm(formula = h2_total ~ Type, family = "binomial", data = .)
# Deviance Residuals:
#   Min        1Q    Median        3Q       Max
# -1.23705  -0.34753  -0.01448   0.34700   1.22204
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.29345    0.01797  -16.33   <2e-16 ***
# TypeMulti    0.55589    0.03461   16.06   <2e-16 ***
  
## Common mQTLs between ARIES and HELIX ####
load("results/eQTLanalysis/comQTLs.Rdata")
mQTLs <- read.table("data/ARIES_mQTLs.tab", header = TRUE, as.is = TRUE)
mQTLsH2 <- read.table("results/ARIES/mqtls.txt", header = TRUE, as.is = TRUE)

## Select pairs detected in previous analysis from ARIES pipeline
## Merge with ARIES original data
ARIESannot <- read.table("data/ariesmqtlsnps.bim", as.is = TRUE)
colnames(ARIESannot) <- c("chr", "SNP", "cm", "pos", "Ref", "Alt")

HELIXannot <- read.table("~/data/WS_HELIX/HELIX_preproc/gwas/Final_data_HRCimp_QC2/HELIX.impQC.rs.bim", as.is = TRUE)
colnames(HELIXannot) <- c("chr", "SNP", "cm", "pos", "Ref", "Alt")

comMQTLs <- mQTLsH2 %>%
  mutate(gene = CpG) %>%
  dplyr::select(SNP, gene, A1, A2, freq, b, se, p, N, r2) %>%
  semi_join(rbind(comCisQTL, comTransQTL), by = c("SNP", "gene")) %>%
  inner_join(left_join(mQTLs, dplyr::select(ARIESannot, SNP, Ref, Alt), by = "SNP"), 
             by = c("SNP", "gene")) %>%
  as_tibble()

## Remove pairs where effect direction do not match
comMQTLs.f <- comMQTLs %>%
  filter(!is.na(Ref)) %>%
  filter(!(sign(b) != sign(beta) & A1 == Ref)) %>%
  filter(!(sign(b) == sign(beta) & A1 == Alt))

comCpGs <- unique(comMQTLs.f$gene)

meQTLTab <- CpGsNum %>%
  mutate(mQTL = CpG %in% comCpGs,
         cisQTL = CpG %in% comCisQTL$gene,
         transQTL = CpG %in% comTransQTL$gene) %>%
  group_by(nCat) %>%
  summarize_if(is.logical, list(sum, mean))

# nCat  mQTL_fn1 cisQTL_fn1 transQTL_fn1 mQTL_fn2 cisQTL_fn2 transQTL_fn2
# <chr>    <int>      <int>        <int>    <dbl>      <dbl>        <dbl>
#  0        26052      24471         1944   0.0742     0.0697      0.00554
#  1         5885       5832          148   0.272      0.270       0.00685
#  2         2427       2408           65   0.335      0.333       0.00898
#  3         1117       1117           24   0.351      0.351       0.00754
#  4+        1141       1131           48   0.356      0.352       0.0150
sum(sigDf$CpG %in% comCpGs)


meQTLTab2 <- CpGsNum %>%
  mutate(mQTLin = CpG %in% comCpGs,
         mQTLOut = !CpG %in% comCpGs) %>%
  group_by(nCat) %>%
  summarize_if(is.logical, sum)

ORs <- sapply(c("1", "2", "3", "4+"), function(x){
  x2 <- data.matrix(subset(meQTLTab2, nCat %in% c("0", x))[, 3:2])
  OR <- x2[1]/x2[2]/x2[3]*x2[4]
  p <- chisq.test(x2)$p.value
  c(OR = OR, p = p)
  
})
meQTL_p <- ORs %>% 
  t() %>%
  data.frame() %>%
  mutate(cat = rownames(.)) %>%
  ggplot(aes(x = cat, y = OR)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "num TCs affected") +
  scale_y_continuous(name = "meQTLs enrichment", trans = "log2") +
  theme_bw() +
  theme(legend.position = "none")

         
         
comMQTLs.f %>%
  select(SNP, gene) %>%
  distinct() %>%
  group_by(gene) %>%
  summarize(CpGs = n()) %>%
  pull(., CpGs) %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.0    13.0    40.0    76.9    98.5  1816.0

comMQTLs.f %>%
  dplyr::select(SNP, gene) %>%
  distinct() %>%
  group_by(gene) %>%
  summarize(CpGs = n()) %>%
  mutate(eQTM = ifelse(gene %in% unique(sigDf$CpG), "eQTM", "non-eQTM")) %>%
  group_by(eQTM) %>%
  summarize(m = median(CpGs),
            l = quantile(CpGs, 0.25),
            h = quantile(CpGs, 0.75))

comMQTLs.f %>%
  dplyr::select(SNP, gene) %>%
  distinct() %>%
  group_by(gene) %>%
  summarize(CpGs = n()) %>%
  mutate(eQTM = ifelse(gene %in% unique(sigDf$CpG), "eQTM", "non-eQTM")) %>%
  glm.nb(CpGs ~ eQTM, .) %>%
  summary()



write.table(meQTLTab[, c(1:2, 5, 3, 6, 4, 7)], file = "paper/meQTL_sum.tab",
            quote = FALSE, row.names = FALSE)

meQTL_p <- CpGsNum %>%
  mutate(mQTL = CpG %in% comCpGs) %>%
  group_by(nCat) %>%
  summarize(prop = mean(mQTL)*100) %>%
  ggplot(aes(x = nCat, y = prop)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_discrete(name = "TCs associated with a CpG") +
  scale_y_continuous(name = "CpGs with meQTL (%)") +
  theme_bw() +
  theme(legend.position = "none")
  

png("paper/eQTMsGenetics.png", width = 3500, height = 2000, res = 300)
plot_grid(plot_grid(h2tot, h2SNP, nrow = 2, labels = "AUTO"), meQTL_p, ncol = 2, labels = c("", "C"))
dev.off()


CpGsNum %>%
  mutate(mQTL = CpG %in% comCpGs) %>%
  glm(mQTL ~ nTC, family = "binomial", .) %>%
  summary()
# Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.391760   0.005913  -404.5   <2e-16 ***
#   nTC          0.518905   0.006034    86.0   <2e-16 ***
  
CpGsNum %>%
  mutate(mQTL = CpG %in% comCpGs) %>%
  filter(nTC > 0) %>%
  glm(mQTL ~ nTC, family = "binomial", .) %>%
  summary()
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.987955   0.017510  -56.42   <2e-16 ***
#   nTC          0.076388   0.006996   10.92   <2e-16 ***
#   



eQTLtab <- CpGsSum %>%
  mutate(mQTL = CpG %in% comCpGs,
         mQTL2 = !mQTL,
         Type = ifelse(Combined != "Non-significant", "Significant", "Non-significant")) %>%
  group_by(Type) %>%
  summarize_if(is.logical, sum) %>%
  ungroup() %>%
  select(-Type) %>%
  data.matrix() 
eQTLtab[2]/eQTLtab[1]/eQTLtab[4]*eQTLtab[3]
# [1] 5.346554


CpGsSum %>%
  filter(Type != "Non-significant") %>%
  mutate(mQTL = CpG %in% comCpGs,
         mQTL2 = !mQTL) %>%
  group_by(Type) %>%
  summarize_if(is.logical, sum) %>%
  ungroup() %>%
  select(-Type) %>%
  data.matrix() %>%
  chisq.test()
# Pearson's Chi-squared test with Yates' continuity correction
# 
# X-squared = 200.52, df = 1, p-value < 2.2e-16
eQTLtab2 <- CpGsSum %>%
  filter(Type != "Non-significant") %>%
  mutate(mQTL = CpG %in% comCpGs,
         mQTL2 = !mQTL) %>%
  group_by(Type) %>%
  summarize_if(is.logical, sum) %>%
  ungroup() %>%
  select(-Type) %>%
  data.matrix() 
eQTLtab2[2]/eQTLtab2[1]/eQTLtab2[4]*eQTLtab2[3]
# [1] 1.39692


# Integrate with gene expression ####
load("results/eQTLanalysis/eQTLs.Rdata")
eQTL <- rbind(gexpme$cis$eqtls, gexpme$trans$eqtls) %>%
  mutate(TC = gene) %>%
  dplyr::select(-gene)

sigDf <- filter(df, sigPair)


## Merge all associations 
## Remove non-coherent associations
mergedDf <- rbind(comCisQTL, comTransQTL) %>%
  mutate(CpG = gene) %>%
  semi_join(comMQTLs.f, by = c("SNP", "gene")) %>%
  inner_join(sigDf, by = "CpG") %>%
  inner_join(eQTL, by = c("snps", "TC")) %>%
  dplyr::select(-sigPair) %>%
  filter(sign(mergedDf$beta.x)*sign(mergedDf$FC) == sign(mergedDf$beta.y))

length(unique(paste(mergedDf$CpG, mergedDf$TC)))
length(unique(paste(mergedDf$CpG, mergedDf$TC)))/nrow(sigDf)
length(unique(mergedDf$CpG))
length(unique(mergedDf$CpG))/length(unique(sigDf$CpG))
length(unique(mergedDf$TC))
length(unique(mergedDf$TC))/length(unique(sigDf$TC))
sum(unique(mergedDf$TC) %in% codingTCs)
sum(unique(mergedDf$TC) %in% codingTCs)/sum(unique(sigDf$TC) %in% codingTCs)


mergedDf <- rbind(comCisQTL, comTransQTL) %>%
  mutate(CpG = gene) %>%
  dplyr::select(-gene, -SNP) %>%
  inner_join(sigDF, by = "CpG") %>%
  inner_join(eQTL, by = c("snps", "TC")) %>%
  dplyr::select(-sigPair)

mergedDf %>%
  select(CpG, TC) %>%
  distinct() %>%
  group_by(TC) %>%
  summarize(CpGs = n()) %>%
  pull(., CpGs) %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   2.000   3.593   4.000  80.000

mergedDf %>%
  select(SNP, TC) %>%
  distinct() %>%
  group_by(TC) %>%
  summarize(CpGs = n()) %>%
  pull(., CpGs) %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.0    19.0    62.0   104.2   145.0  1114.0



mergedDf %>%
  select(SNP, CpG) %>%
  distinct() %>%
  group_by(CpG) %>%
  summarize(CpGs = n()) %>%
  pull(., CpGs) %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00   17.00   51.00   88.27  122.00 1286.00

mergedDf %>%
  select(TC, CpG) %>%
  distinct() %>%
  group_by(CpG) %>%
  summarize(CpGs = n()) %>%
  pull(., CpGs) %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00    1.00    1.00    1.87    2.00   21.00


## Check replication in GTEx
### Download data (https://gtexportal.org/home/datasets) - 28/11/2019
tar --extract --file=data/GTEx_Analysis_v8_eQTL.tar GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz
mv GTEx_Analysis_v8_eQTL/ data/
  

### Example
# SNP: rs11585123
# CpG: cg15580684
# Gene: AJAP1 (TC01000080.hg.1)

plink <- read.plink("/home/isglobal.lan/cruiz/data/WS_HELIX/HELIX_preproc/gwas/Final_data_HRCimp_QC2/HELIX.impQC.rs")
load("results/preprocessFiles/Methylation_GRSet.RData")
colnames(gset) <- gset$HelixID

load("results/preprocessFiles/Expression_SE_residuals.RData")
colnames(se) <- se$HelixID

## Common samples
comSamps <- intersect(intersect(rownames(plink$geno), colnames(gset)), colnames(se))

alleles <- plink$map["rs11585123", c("allele.1", "allele.2")]
alleles <- c(paste0(alleles[1], alleles[1]),
             paste0(alleles[1], alleles[2]),
             paste0(alleles[2], alleles[2]))

dat <- data.frame(geno = as.numeric(plink$geno[comSamps, "rs11585123"]),
                  methy = as.numeric(getBeta(gset["cg15580684", comSamps])),
                  gexp = as.numeric(assay(se["TC01000080.hg.1", comSamps])))
dat$geno <- alleles[dat$geno]

sm <- ggplot(dat, aes(x = factor(geno), y = methy)) + 
  geom_boxplot() +
  scale_x_discrete(name = "rs11585123") +
  scale_y_continuous(name = "cg15580684") +
  theme_bw()

sg <- ggplot(dat, aes(x = factor(geno), y = gexp)) + 
  geom_boxplot() +
  scale_x_discrete(name = "rs11585123") +
  scale_y_continuous(name = "TC01000080.hg.1") +
  theme_bw()

  
me <- ggplot(dat, aes(x = methy, y = gexp)) + 
  geom_point() +
  scale_x_continuous(name = "cg15580684") +
  scale_y_continuous(name = "TC01000080.hg.1") +
  geom_smooth(method = "lm") +
  theme_bw()
png("paper/eQTMstrio.png", width = 3500, height = 2000, res = 300)
plot_grid(plot_grid(sm, sg, nrow = 1), me, nrow = 2)
dev.off()

## Run enrichment analysis genes in trios ####
## Copy functions from eQTM_interpretation.R
snpGenes <- df %>%
  dplyr::select(TC) %>%
  distinct() %>%
  mutate(sig = factor(ifelse(TC %in% unique(mergedDf$TC), 1, 0)))
snpGos <- computeGOs(snpGenes)
save(snpGos, file = "paper/snpGOs.Rdata")