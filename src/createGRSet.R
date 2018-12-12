###############################################################################
#' Create GenomicRatioSet for analysis (methylation)
###############################################################################

library('minfi')

# Loading GRSet
load("/Lacie_CRW10023/HELIX_preproc/methylation/Final_data/methylome_subcohort_ComBatSlide_6cells_v3.Rdata")
gset <- methylome_subcohort_ComBatSlide_6cells
pheno <- pData(gset)

# ---- Pheno
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
dim(merg) # 933 11

# Caucasians only filtering
pheno <- pheno[pheno$h_ethnicity_cauc == "yes",]
dim(pheno)
table(pheno$h_ethnicity_cauc)

# Sample IDs filtering on GRset
filter <- pheno$SampleID %in% merg$SampleID
pheno <- pheno[filter,]
dim(pheno)

ids <- rownames(pheno)
save(ids, file = "ids.RData")

# Adding the IDs to the GRset
gset <- gset[,ids]

# ---- Annotation
# Loading

rd <- getAnnotation(gset)
rd <- rd[,-c(5:14, 16, 20:21)]

# Adding phantom summary
rd$Phantom_S <- ifelse(grepl("low", rd$Phantom), "low",
                         ifelse(grepl("high", rd$Phantom),"high", ""))

# Adding group columns
rd$TSS200 <- ifelse(grepl("TSS200", rd$UCSC_RefGene_Group), T, F)
rd$TSS1500 <- ifelse(grepl("TSS1500", rd$UCSC_RefGene_Group), T, F)
rd$UTR5 <- ifelse(grepl("5'UTR", rd$UCSC_RefGene_Group), T, F)
rd$FirstExon <- ifelse(grepl("1stExon", rd$UCSC_RefGene_Group), T, F)
rd$Body <- ifelse(grepl("Body", rd$UCSC_RefGene_Group), T, F)
rd$UTR3 <- ifelse(grepl("3'UTR", rd$UCSC_RefGene_Group), T, F)

# Adding dhs and crom15
crom15 <- read.csv("crom15.csv", header = T, row.names = 1)
crom15 <- crom15[,c(1,7:21)]
dhs <- read.csv("dhs.csv", header = T, row.names = 1)
dhs <- dhs[,c(1,7:19)]

states <- merge(crom15, dhs, by = "HT12v4.ArrayAddress", sort = F)
rd <- merge(rd, states, by.x = "row.names", by.y = "HT12v4.ArrayAddress", sort = F)

a <- strsplit(rd$UCSC_RefGene_Name, ";")
a[a == "character(0)"] <- ""
rd$UCSC_RefGene_Name <- a


# Applying new rowData
rowData(gset) <- rd

# ---- Saving all
# save(gset, file = "gset_sm_1.RData")
# load("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/gset_sm_1.RData")



