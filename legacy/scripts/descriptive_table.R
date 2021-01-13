###############################################################
# DESCRIPTIVES TABLE
###############################################################

# Transcriptome ID list loading
tran <- read.table("/Lacie_CRW10023/HELIX_preproc/gene_expression/Final_data/ID_list_transcriptome_v2.txt", header = T)
head(tran)

# Period 1B filtering
tran <- tran[tran$Period != "1B",]
nrow(tran)
table(tran$Period)

# Methylome ID list loading
met <- read.table("/Lacie_CRW10023/HELIX_preproc/methylation/Final_data/ID_list_methylome_v3.txt", header = T)
head(met)

# Period IB filtering
met <- met[met$Period != "1B",]
nrow(met)
table(met$Period)

# Transcriptome and methylome IDs merging
merg <- merge(tran, met, by = "SampleID")
dim(merg)

# GenomicRatioSet loading
library("minfi")
load("/Lacie_CRW10023/HELIX_preproc/methylation/Final_data/methylome_subcohort_ComBatSlide_6cells_v3.Rdata")
gset <- methylome_subcohort_ComBatSlide_6cells
pheno <- pData(gset)
head(pheno)
dim(pheno)

# Caucasians only filtering
pheno <- pheno[pheno$h_ethnicity_cauc == "yes",]
dim(pheno)
table(pheno$h_ethnicity_cauc)

# Sample IDs filtering on GRset
filter <- pheno$SampleID %in% merg$SampleID
pheno <- pheno[filter,]
dim(pheno)

###############################################################

names(pheno)

# Cohort
table(pheno$cohort)
prop.table(table(pheno$cohort))*100

# Age
mean(pheno$age_sample_years)
sd(pheno$age_sample_years)

# Sex
table(pheno$e3_sex)
prop.table(table(pheno$e3_sex))*100

# Ethnicity omitted (all caucasians)

# Child BMI
mean(pheno$hs_zbmi_theano)
sd(pheno$hs_zbmi_theano)

# NK prop
mean(pheno$NK_6)
sd(pheno$NK_6)

# Bcell prop
mean(pheno$Bcell_6)
sd(pheno$Bcell_6)

# CD4T prop
mean(pheno$CD4T_6)
sd(pheno$CD4T_6)

# CD8T prop
mean(pheno$CD8T_6)
sd(pheno$CD8T_6)

# Gran prop
mean(pheno$Gran_6)
sd(pheno$Gran_6)

# Mono prop
mean(pheno$Mono_6)
sd(pheno$Mono_6)

# Time from blood extraction to last meal
mean(pheno$hs_dift_mealblood_imp)
sd(pheno$hs_dift_mealblood_imp)

# Time of the day when blood sample was collected
mean(pheno$blood_sam4)
sd(pheno$blood_sam4)
