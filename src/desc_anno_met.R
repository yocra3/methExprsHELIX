###############################################################
# ANNOTATION DESCRIPTIVES (METHYLATION)
###############################################################

# GenomicRatioSet loading
library("minfi")
load("/Lacie_CRW10023/HELIX_preproc/methylation/Final_data/methylome_subcohort_ComBatSlide_6cells_v3.Rdata")
gset <- methylome_subcohort_ComBatSlide_6cells
anno <- getAnnotation(gset)
head(anno)
dim(anno)

# Group table
group <- strsplit(anno$UCSC_RefGene_Group, ";")
group <- sapply(group, head, n=1)
group <- unlist(lapply(group,function(x) if(identical(x,character(0))) ' ' else x))
table(group)

# Group table (totals)
# group.s <- paste(anno$UCSC_RefGene_Group, collapse = ";")
# group.v <- strsplit(group.s, ";")[[1]]
# table(group.v)

#
prom <- anno[grep("TSS1500|TSS200|5'UTR|1stExon", anno$UCSC_RefGene_Group),]
dim(prom)

# Relation to Island
table(anno$Relation_to_Island)

# Intergenic
inter <- anno[anno$UCSC_RefGene_Name == "",]
nrow(inter)
#Genic
genic <- anno[anno$UCSC_RefGene_Name != "",]
nrow(genic)

# CpGs & amount of genes
genic <- anno[anno$UCSC_RefGene_Name != "",]
cpgenes <- genic$UCSC_RefGene_Name
cpgenes2 <- strsplit(cpgenes, ";")
counts <- lapply(cpgenes2, length)
counts <- unlist(counts)
table(counts)
# write.csv(counts, "counts2.csv")

# Chr and amount of CpGs
chr <- as.data.frame(table(anno$chr))
number <- gsub("chr", "", chr$Var1)
number <- as.numeric(number)
chr <- cbind(number, chr)
chr <- chr[order(chr$number),]
# write.csv(chr, "chr.csv")

# Summarizing Phantom
anno$Phantom_S <- ifelse(grepl("low", anno$Phantom), "low", ifelse(grepl("high", anno$Phantom),"high", ""))
table(anno$Phantom_S)

# Deleting some columns
anno <- anno[ ,-which(names(anno) %in% c("AddressA", "AddressB", "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color", "Forward_Sequence", "SourceSeq", "Random_Loci", "Methyl27_Loci"))]

#################################

# DMR
table(anno$DMR)

# Group columns
anno$TSS200 <- ifelse(grepl("TSS200", anno$UCSC_RefGene_Group), T, F)
anno$TSS1500 <- ifelse(grepl("TSS1500", anno$UCSC_RefGene_Group), T, F)
anno$UTR5 <- ifelse(grepl("5'UTR", anno$UCSC_RefGene_Group), T, F)
anno$FirstExon <- ifelse(grepl("1stExon", anno$UCSC_RefGene_Group), T, F)
anno$Body <- ifelse(grepl("Body", anno$UCSC_RefGene_Group), T, F)
anno$UTR3 <- ifelse(grepl("3'UTR", anno$UCSC_RefGene_Group), T, F)

groupdf <- cbind("TSS200" = anno$TSS200,"TSS1500" = anno$TSS1500,
                 "5'UTR" = anno$UTR5, "1stExon" = anno$FirstExon,
                 "Body" = anno$Body, "3'UTR" = anno$UTR3)

rownames(groupdf) <- anno$Name
write.csv(groupdf, "groups.csv")

# CpGs & Genes (without duplicates)
genic <- anno[anno$UCSC_RefGene_Name != "",]
cpgenes <- genic$UCSC_RefGene_Name
cpgenes2 <- strsplit(cpgenes, ";")
cpgenes3 <- lapply(cpgenes2, unique)
counts <- lapply(cpgenes3, length)
counts <- unlist(counts)
counts <- as.data.frame(table(counts))
colnames(counts) <- c("N_Genes", "CpG_freq")
# write.csv(counts, "cpg_genes.csv")

# Genes (nom) - CpGs
gene.s <- paste(anno$UCSC_RefGene_Name, collapse = ";")
gene.v <- strsplit(gene.s, ";")[[1]]
genes <- unique(gene.v)

cpgcount <- function(x, arg1) {
  sum(grepl(x, arg1))
}

cpgc <- lapply(genes, cpgcount, arg1 = anno$UCSC_RefGene_Name)
cpgc2 <- unlist(cpgc)
genes_cpg <- data.frame("Genes" = genes, "CpG_count" = cpgc2)
write.csv(genes_cpg, "genes_cpgc.csv")

######################