library(biomaRt)

#        BroadTUCP          ENSEMBL          GenBank       GenBankHTC
#             7203            90058            51431             5162
# Havanatranscript         lncRNADb          NONCODE           RefSeq
#            92250              353            47175            37802
#      RinnlincRNA        UCSCGenes
#              714            75072

############################################################
# Look for TSS in BioMart

load("/Lacie_CRW10023/HELIX_preproc/gene_expression/annotation/TC_mrna_ID_DB.RData")
TC_mrna_ID_DB$TSS_Affy <- ifelse(TC_mrna_ID_DB$strand == "+",
                                 TC_mrna_ID_DB$start, TC_mrna_ID_DB$stop)

# refseq. ok
# UCSCGenes. ok (link to anther db refseq) -> xxxx.
# ENSEMBL. ok (link to anther db refseq) -> chr, start, end, gene symbol ... quants???
# Havanatranscript. ... biomart
# 
# GenBank
# GenBankHTC
# 
# BroadTUCP
# lncRNADb          
# RinnlincRNA        
# noncode

# List of mrna IDs
mrnas <- TC_mrna_ID_DB$mrna_ID
any(duplicated(mrnas))

# Take into acount that there are duplicates (?)
c <- TC_mrna_ID_DB[!duplicated(mrnas),]
table(c$mrna_DB)

#          BroadTUCP          ENSEMBL          GenBank       GenBankHTC
#               7203            89887            47467             5063
#   Havanatranscript         lncRNADb          NONCODE           RefSeq
#              92014               40            31831            36824
#        RinnlincRNA        UCSCGenes
#                714            74907

# # # (!) NONCODE: long-non-coding ids from different databases 

library(biomaRt)
grch37 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",
                  path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")

## ENSEMBL
d1 <- getBM(attributes = c('ensembl_transcript_id','transcription_start_site',
                            'hgnc_symbol', 'chromosome_name', 'start_position',
                            'end_position', 'strand'),
            filters = 'ensembl_transcript_id', 
            values = mrnas, 
            mart = grch37)

dim(d1) # 83627     7
# From 89887 to 83627
d1$DB <- "ENSEMBL"

## Havana Transcript
d2 <- getBM(attributes = c('ottt','transcription_start_site',
                            'hgnc_symbol', 'chromosome_name', 'start_position',
                            'end_position', 'strand'),
            filters = 'ottt', 
            values = mrnas, 
            mart = grch37)

dim(d2) # 89895     7
# From 92014 to 89895
d2$DB <- "Havanatranscript"

## UCSC (via web)
noncode_id <- c$mrna_ID[which(c$mrna_DB == "NONCODE")]
uc_id1 <- noncode_id[grep("uc", noncode_id)]
uc_id2 <- c$mrna_ID[which(c$mrna_DB == "UCSCGenes")]
uc_id <- c(uc_id1, uc_id2)
length(uc_id) # 77247

write(uc_id, file = "uc_id.txt", sep = "\n")
# File "uc_id.txt" was uploaded to UCSC Table Browser

# Web warning:
#   Note: 2543 of the 77247 given identifiers have no match in table
#   knownGene, field name or in alias table kgAlias, field alias.

# Downloaded file:
#   "UCSCGenes_TSS_hg19.txt"

d3 <- read.delim("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/UCSCGenes_TSS_hg19.txt")
dim(d3) # 76157     9

table(d3$X.hg19.knownGene.name %in% uc_id)
# FALSE  TRUE
#  7145 69012

# Not all the uc_id match to the downloaded table
# Fixing:
d3 <- d3[d3$X.hg19.knownGene.name %in% uc_id,]
dim(d3) # 69012     9

## RefSeq
d4 <- getBM(attributes = c('refseq_mrna','transcription_start_site',
                           'hgnc_symbol', 'chromosome_name', 'start_position',
                           'end_position', 'strand'),
            filters = 'refseq_mrna', 
            values = unique(mrnas), 
            mart = grch37)

dim(d4) # 32111     7
d4$DB <- "Refseq_mrna"

d5 <- getBM(attributes = c('refseq_ncrna','transcription_start_site',
                           'hgnc_symbol', 'chromosome_name', 'start_position',
                           'end_position', 'strand'),
            filters = 'refseq_ncrna', 
            values = mrnas, 
            mart = grch37)

dim(d5) # 3544    7
d5$DB <- "Refseq_ncrna"

# From 36824 to 35655 (32111 + 3544)
# Note: Predicted mrna id and ncrna id attributes returned 0 rows

## Organizing obtained data

#UCSC
d3 <- d3[,c(1,4,9,2,6,7,3)]
dim(d3) # 76157     7

# save(d3, file = "d3.RData")


colnames(d1)[1] <- "ID"
colnames(d2)[1] <- "ID"
colnames(d4)[1] <- "ID"
colnames(d5)[1] <- "ID"
d1245 <- rbind(d1, d2, d4, d5)
d1245$strand <- ifelse(d1245$strand == "1", "+", "-")
dim(1245) # 209177      7

# save(d1245, file = "d1245.RData")

# ----
load("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/TC_mrna_ID_DB.RData")
load("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/d3.RData")
load("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/d1245.RData")

head(d3)
d3$DB <- "UCSC"
head(d1245)
d3[, c(1,2,3,8,9)] <- lapply(d3[, c(1,2,3,8,9)], as.character)

table(d3$hg19.knownGene.chrom)
table(d1245$chromosome_name)

hf <- grep("hap", d3$hg19.knownGene.chrom)
d3 <- d3[-(hf),]
hf <- grep("random", d3$hg19.knownGene.chrom)
d3 <- d3[-(hf),]

hf <- grep("H", d1245$chromosome_name)
d1245 <- d1245[-(hf),]

table(d3$hg19.knownGene.chrom)
table(d1245$chromosome_name)

s <- sapply(d3$hg19.knownGene.chrom, function(x) gsub("chr", "", x))
d3$hg19.knownGene.chrom <- s
table(d3$hg19.knownGene.chrom)

prep1 <- d3[,c(1,2,3,4,5,9)]
colnames(prep1) <- c("ID", "chromDB", "strandDB", "startDB", "endDB",
                     "GeneSymbolDB")
prep1$TSSDB <- ifelse(prep1$strand == "+", prep1$startDB, prep1$endDB)
prep1$DB <- "UCSC"

prep2 <- d1245[,c(1,4,7,5,6,3,2,8)]
colnames(prep2) <- c("ID", "chromDB", "strandDB", "startDB", "endDB",
                     "GeneSymbolDB", "TSSDB", "DB")

head(prep1)
head(prep2)
dfinal <- rbind(prep2, prep1)
# save(dfinal, file = "dfinal.RData")
# load("dfinal.RData")

ndf <- merge(TC_mrna_ID_DB, dfinal, by.x = "mrna_ID", by.y = "ID")
ndf$tssdif <- ndf$TSSDB - ndf$TSS_Affy
ndf$tssdif <- ifelse(ndf$strand == "+", ndf$tssdif, -ndf$tssdif)
dim(ndf) # 286109     23

ndf$tssdif2 <- NULL
ndf$tssdif2[ndf$tssdif == 0] <- "0"
ndf$tssdif2[ndf$tssdif < 0] <- "-1:-100"
ndf$tssdif2[ndf$tssdif < -100] <- "-100:-1K"
ndf$tssdif2[ndf$tssdif < -1000] <- "-1K:-10K"
ndf$tssdif2[ndf$tssdif < -10000] <- "-10K:-100K"
ndf$tssdif2[ndf$tssdif < -100000] <- "<-100K"
ndf$tssdif2[ndf$tssdif > 0] <- "1:100"
ndf$tssdif2[ndf$tssdif > 100] <- "100:1K"
ndf$tssdif2[ndf$tssdif > 1000] <- "1K:10K"
ndf$tssdif2[ndf$tssdif > 10000] <- "10K:100K"
ndf$tssdif2[ndf$tssdif > 100000] <- ">100K"

levs <- c("<-100K", "-10K:-100K", "-1K:-10K", "-100:-1K", "-1:-100", "0",
          "1:100", "100:1K", "1K:10K", "10K:100K", ">100K")
ndf$tssdif2 <- as.factor(ndf$tssdif2)
ndf$tssdif2 <- factor(ndf$tssdif2, levels = levs)

ndf$Multiple <- ifelse(ndf$GeneSymbol_Affy == "", "None",
                       ifelse(grepl(";", ndf$GeneSymbol_Affy), "Multiple", "One"))
ndf$Multiple <- as.character(ndf$Multiple)

tcc <- as.data.frame(table(ndf$TC))
dim(tcc) # 45004     2

ndf <- merge(ndf, tcc, by.x = "TC", by.y = "Var1")
dim(ndf) # 286109     26
# save(ndf, file = "ndf.RData")

#####################################################################
# PLOTS

# load("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/ndf.RData")

ndf <- ndf[ndf$strand == ndf$strandDB,]
# 285178 rows

library(ggplot2)
ndf$Freq2 <- NULL
ndf$Freq2[ndf$Freq == 0] <- "0"
ndf$Freq2[ndf$Freq == 1] <- "1"
ndf$Freq2[ndf$Freq > 1] <- "2"
ndf$Freq2[ndf$Freq > 2] <- "3-5"
ndf$Freq2[ndf$Freq > 5] <- "6-10"
ndf$Freq2[ndf$Freq > 10] <- "11-20"
ndf$Freq2[ndf$Freq > 20] <- "21-40"
ndf$Freq2[ndf$Freq > 40] <- "41-60"
ndf$Freq2[ndf$Freq > 60] <- "61-80"
ndf$Freq2[ndf$Freq > 80] <- ">80"

levs2 <- c("0", "1", "2", "3-5", "6-10", "11-20", "21-40", "41-60",
           "61-80", ">80")
ndf$Freq2 <- as.factor(ndf$Freq2)
ndf$Freq2 <- factor(ndf$Freq2, levels = levs2)

##################
levs <- c("<-100K", "-10K:-100K", "-1K:-10K", "-100:-1K", "-1:-100", "0",
          "1:100", "100:1K", "1K:10K", "10K:100K", ">100K")
ndf$tssdif2 <- factor(ndf$tssdif2, levels = rev(levs))

p1 <- ndf[ndf$strand == "+",]
p1 <- as.data.frame(prop.table(table(p1$tssdif2)))
p1 <- round(p1[,2] * 100, 2)

p2 <- ndf[ndf$strand == "-",]
p2 <- as.data.frame(prop.table(table(p2$tssdif2)))
p2 <- round(p2[,2] * 100, 2)

p3 <- ndf[ndf$phase == "phase1",]
p3 <- as.data.frame(prop.table(table(p3$tssdif2)))
p3 <- round(p3[,2] * 100, 2)

p4 <- ndf[ndf$phase == "phase2",]
p4 <- as.data.frame(prop.table(table(p4$tssdif2)))
p4 <- round(p4[,2] * 100, 2)

p5 <- ndf[ndf$Multiple == "Multiple",]
p5 <- as.data.frame(prop.table(table(p5$tssdif2)))
p5 <- round(p5[,2] * 100, 2)

p6 <- ndf[ndf$Multiple == "None",]
p6 <- as.data.frame(prop.table(table(p6$tssdif2)))
p6 <- round(p6[,2] * 100, 2)

p7 <- ndf[ndf$Multiple == "One",]
p7 <- as.data.frame(prop.table(table(p7$tssdif2)))
p7 <- round(p7[,2] * 100, 2)

p8 <- ndf[ndf$DB == "ENSEMBL",]
p8 <- as.data.frame(prop.table(table(p8$tssdif2)))
p8 <- round(p8[,2] * 100, 2)

p9 <- ndf[ndf$DB == "Havanatranscript",]
p9 <- as.data.frame(prop.table(table(p9$tssdif2)))
p9 <- round(p9[,2] * 100, 2)

p10 <- ndf[ndf$DB == "Refseq_mrna",]
p10 <- as.data.frame(prop.table(table(p10$tssdif2)))
p10 <- round(p10[,2] * 100, 2)

p11 <- ndf[ndf$DB == "Refseq_ncrna",]
p11 <- as.data.frame(prop.table(table(p11$tssdif2)))
p11 <- round(p11[,2] * 100, 2)

p12 <- ndf[ndf$DB == "UCSC",]
p12 <- as.data.frame(prop.table(table(p12$tssdif2)))
p12 <- round(p12[,2] * 100, 2)

p13 <- ndf[ndf$Freq2 == "1",]
p13 <- as.data.frame(prop.table(table(p13$tssdif2)))
p13 <- round(p13[,2] * 100, 2)

p14 <- ndf[ndf$Freq2 == "2",]
p14 <- as.data.frame(prop.table(table(p14$tssdif2)))
p14 <- round(p14[,2] * 100, 2)

p15 <- ndf[ndf$Freq2 == "3-5",]
p15 <- as.data.frame(prop.table(table(p15$tssdif2)))
p15 <- round(p15[,2] * 100, 2)

p16 <- ndf[ndf$Freq2 == "6-10",]
p16 <- as.data.frame(prop.table(table(p16$tssdif2)))
p16 <- round(p16[,2] * 100, 2)

p17 <- ndf[ndf$Freq2 == "11-20",]
p17 <- as.data.frame(prop.table(table(p17$tssdif2)))
p17 <- round(p17[,2] * 100, 2)

p18 <- ndf[ndf$Freq2 == "21-40",]
p18 <- as.data.frame(prop.table(table(p18$tssdif2)))
p18 <- round(p18[,2] * 100, 2)

p19 <- ndf[ndf$Freq2 == "41-60",]
p19 <- as.data.frame(prop.table(table(p19$tssdif2)))
p19 <- round(p19[,2] * 100, 2)

p20 <- ndf[ndf$Freq2 == "61-80",]
p20 <- as.data.frame(prop.table(table(p20$tssdif2)))
p20 <- round(p20[,2] * 100, 2)

p21 <- ndf[ndf$Freq2 == ">80",]
p21 <- as.data.frame(prop.table(table(p21$tssdif2)))
p21 <- round(p21[,2] * 100, 2)

ndfp <- data.frame("tssdif2" = levels(ndf$tssdif2),
                   "strand_plus" = p1, "strand_minus" = p2,
                   "phase1" = p3, "phase2" = p4,
                   "Multiple_gene" = p5, "None_genes" = p6,
                   "One_gene" = p7, "ENSEMBL" = p8, "Havanatranscript" = p9,
                   "Refseq_mrna" = p10, "Refseq_ncrna" = p11, "UCSC" = p12,
                   "x1_mrna" = p13, "x2_mrnas" = p14, "x3_5_mrnas" = p15,
                   "x6_10_mrnas" = p16, "x11_20_mrnas" = p17,
                   "x21_40_mrnas" = p18, "x41_60_mrnas" = p19,
                   "x61_80_mrnas" = p20, "more_than_80_mrnas" = p21)

levs <- c("<-100K", "-10K:-100K", "-1K:-10K", "-100:-1K", "-1:-100", "0",
          "1:100", "100:1K", "1K:10K", "10K:100K", ">100K")

ndfp$tssdif2 <- factor(ndfp$tssdif2, levels = rev(levs))
prova <- melt(ndfp, id = "tssdif2")

pdf("ndf.pdf", height = 5, width = 12)

ggplot(data = prova, aes(x = variable, y = value, fill = tssdif2)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_text(aes(label = paste(value,"%",sep="")),size = 2, vjust = 1.5, position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

pdf("ndf2.pdf", height = 6, width = 12)

# Without logs
ggplot(ndf, aes(tssdif)) +
  geom_density() +
  ggtitle("All") +
  xlim(-100, 500)

ggplot(ndf, aes(tssdif, color = strand)) +
  geom_density() +
  ggtitle("Strand") +
  xlim(-100, 500)

ggplot(ndf, aes(tssdif, color = phase)) +
  geom_density() +
  ggtitle("Phase") +
  xlim(-100, 500)

ggplot(ndf, aes(tssdif, color = geneNA)) +
  geom_density() +
  ggtitle("Gene NA") +
  xlim(-100, 500)

ggplot(ndf, aes(tssdif, color = Multiple)) +
  geom_density() +
  ggtitle("Multiple") +
  xlim(-100, 500)

ggplot(ndf, aes(tssdif, color = DB)) +
  geom_density() +
  ggtitle("Database") +
  xlim(-100, 500)

ggplot(ndf, aes(tssdif, color = Freq2)) +
  geom_density() +
  ggtitle("mRNA count") +
  xlim(-100, 500)

# With logs
ggplot(ndf, aes(log10(tssdif))) +
  geom_density(alpha = 0.5) +
  ggtitle("All")

ggplot(ndf, aes(log10(tssdif), color = strand)) +
  geom_density(alpha = 0.5) +
  ggtitle("Strand")

ggplot(ndf, aes(log10(tssdif), color = phase)) +
  geom_density(alpha = 0.5) +
  ggtitle("Phase")

ggplot(ndf, aes(log10(tssdif), color = geneNA)) +
  geom_density(alpha = 0.5) +
  ggtitle("Gene NA")

ggplot(ndf, aes(log10(tssdif), color = Multiple)) +
  geom_density(alpha = 0.5) +
  ggtitle("Multiple")

ggplot(ndf, aes(log10(tssdif), color = DB)) +
  geom_density(alpha = 0.5) +
  ggtitle("Database")

ggplot(ndf, aes(log10(tssdif), color = Freq2)) +
  geom_density(alpha = 0.5) +
  ggtitle("mRNA count")

dev.off()

#####################################################################
# COMPLEXING

load("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/ndf.RData")
load("/Lacie_CRW10023/HELIX_analyses/expr_met_SM/TC_mrna_ID_DB.RData")
load("/Lacie_CRW10023/HELIX_preproc/gene_expression/annotation/HTAv2_0.na36.hg19.transcript_filtr.RData")
dim(ndf) #286109  26
dim(TC_mrna_ID_DB) # 407337     15
anno <- HTAv2_0.na36.hg19.transcript_filtr
dim(anno) # 64828    15

tcc <- as.data.frame(table(ndf$TC))
tcs1 <- as.character(tcc$Var1)

tcc <- as.data.frame(table(TC_mrna_ID_DB$TC))
tcs2 <- as.character(tcc$Var1)

cx1 <- function(x){
  d <- ndf[,"GeneSymbolDB"]
  vec <- d[which(ndf$TC == x)]
}

cx2 <- function(x, coln){
  filt <- which(TC_mrna_ID_DB$TC == x)
  vec <- TC_mrna_ID_DB[filt, coln]
}

GeneSymbolDB <- lapply(tcs1, cx1)
GeneSymbolDB2 <- lapply(GeneSymbolDB, unique)
GeneSymbolDB2 <- lapply(GeneSymbolDB2, function(x){x[x!=""]})
GeneSymbolDB2[GeneSymbolDB2 == "character(0)"] <- ""
mrna_ID <- lapply(tcs2, cx2, coln = "mrna_ID")
mrna_DB <- lapply(tcs2, cx2, coln = "mrna_DB")


temp1 <- data.frame("transcript_cluster_id" = tcs1)
temp1$GeneSymbolDB <- GeneSymbolDB
temp1$GeneSymbolDB2 <- GeneSymbolDB2

temp2 <- data.frame("transcript_cluster_id" = tcs2)
temp2$mrna_ID <- mrna_ID
temp2$mrna_DB <- mrna_DB

anno1 <- merge(anno, temp1, by = "transcript_cluster_id", all.x = TRUE)
anno1 <- merge(anno1, temp2, by = "transcript_cluster_id")

cm <- function(x){
  if (!is.na(x)){length(x)}
  else {0}
}

mrna_N <- lapply(anno1$mrna_ID, cm)
anno1$mrna_N <- mrna_N

x <- strsplit(anno1$GeneSymbol_Affy, ";")
x[x == "character(0)"] <- ""
anno1$GeneSymbol_Affy <- x

HTAv2_0.na36.hg19.transcript_filtr_2 <- anno1
save(HTAv2_0.na36.hg19.transcript_filtr_2, file = "HTAv2_0.na36.hg19.transcript_filtr_2.RData")

#####################################################################
# SUMMARIZED EXPERIMENT

load("HTAv2_0.na36.hg19.transcript_filtr_2.RData")
anno1 <- HTAv2_0.na36.hg19.transcript_filtr_2
rownames(anno1) <- anno1$transcript_cluster_id
dim(anno1) # 64828    20

# Loading transcriptome
library('Biobase')
#load("/Lacie_CRW10023/HELIX_preproc/gene_expression/Final_data/transcriptome_subcohort_f1_v2.Rdata")
load("/scratch/smari/transcriptome_subcohort_f1_v2.Rdata")
eset <- transcriptome_subcohort_f1
feat <- featureData(eset)
eanno <- pData(feat)
names(eanno)
dim(eanno) # 59439    21


# Filtering undetermined chromosome positions & mitochondrial
eanno$chromosome <- as.character(eanno$chromosome)
eanno <- eanno[-grep("gl", eanno$chromosome),]
eanno <- eanno[-grep("chrM", eanno$chromosome),]
dim(eanno) # 59343    21

# Intersecting TCs
TCs <- intersect(eanno$transcript_cluster_id, anno1$transcript_cluster_id)
anno1 <- anno1[TCs,]
dim(anno1) # 58501    20

# Filter TCs
eset <- eset[TCs,]
dim(eset)
# Features  Samples
# 58501     1161

# Adding anno
pData(feat) <- anno1
featureData(eset) <- feat

# Filter IDs (from methylomes)
load("/scratch/smari/ids.RData")
eset <- eset[,ids]
dim(eset)
# Features  Samples
# 58501      832

se <- makeSummarizedExperimentFromExpressionSet(eset)
validObject(se) # [1] TRUE
save(se, file = "summ_exp_sm_1.RData")

# Download chr 22
gr <- rowRanges(se)
se22 <- se[seqnames(gr) == 'chr22', ]
save(se22, file = "se22.RData")

#####################################################################
# FINDOVERLAPS

# Residuals made in code_residuals2.R

library(minfi)
library(ggplot2)

# load("summ_exp_sm_1.RData")
load("summ_exp_sm_res1.RData")
# se
# load("gset_sm_1.RData")
load("gset_sm_res1.RData")
# gset

# Changing colname 4 to prevent findOverlaps error
rdgset <- rowData(gset_res)
colnames(rdgset)[4] <- "str"
rowData(gset_res) <- rdgset

fo <- findOverlaps(rowRanges(gset_res) + 5e5, rowRanges(se_res))
# 13,973,780 hits

froms <- from(fo)
tos <- to(fo) 
descfrom <- as.data.frame(table(froms))
descto <- as.data.frame(table(tos))

cpg <- rownames(gset_res)
tc <- rownames(se_res)

froms2 <- cpg[froms]
tos2 <- tc[tos]


pdf("dens_fo.pdf", width = 8, height = 6)

ggplot(descfrom, aes(x = Freq)) +
  geom_histogram(color = "steelblue", fill = "white") +
  scale_x_continuous(breaks = c(0, 50, 100 ,150, 200)) +
  theme_classic()

ggplot(descto, aes(x = Freq)) +
  geom_histogram(color = "seagreen", fill = "white") +
  theme_classic()

dev.off()

pairs <- data.frame("CpG" = froms2, "TC" = tos2)
save(as.matrix(pairs), file = "pairs.RData")

gr_met <- as.data.frame(rowRanges(gset_res))
gr_se <- as.data.frame(rowRanges(se_res))

pairs_data <- merge(pairs, gr_se[,c(1:5,16)], by.x = "TC", by.y = "row.names")
pairs_data <- merge(pairs_data, gr_met[,1:5], by.x = "CpG", by.y = "row.names")
pairs_data$dif <- pairs_data$start.y - pairs_data$TSS_Affy
pairs_data$seqnames.x <- as.character(pairs_data$seqnames.x)
pairs_data$seqnames.y <- as.character(pairs_data$seqnames.y)

summary(pairs_data$dif)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2935745.0  -241866.0       60.0      193.8   241444.0  2798370.0
table(pairs_data$seqnames.x == pairs_data$seqnames.y)
# 
# TRUE
# 13973780

save(pairs_data, file = "pairs_data.RData")

# Only chromosome 22
load('pairs_data.RData')
pairs <- pairs_data[pairs_data$seqnames.x == 'chr22',]
pairs <- as.matrix(pairs[,1:2])
save(pairs, file = "pairs22.RData")

eassay <- assays(se_res)$exprs
massay <- assays(gset_res)$Beta
save(eassay, file = "eassay_res1.RData")
save(massay, file = "massay_res1.RData")

# linear_m.R

# dataframe amb noms, no numeros (ja guardat tot)
# pasar com arguments maassay i eassay
# escrigui un fitxer amb la matriu sense colnames
# en lloc de load se i gset, nomes les matrius
# mirar quan ocupa findoverlaps amb noms posats
# provar correr script amb cromosoma 22
#   Rscript fitxer arguments
#   Mesurar temps
#   provar amb chr 1
#   paralel mcapply 8 cores