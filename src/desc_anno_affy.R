###############################################################
# ANNOTATION DESCRIPTIVES (EXPRESSION)
###############################################################

library('Biobase')

load("/Lacie_CRW10023/HELIX_preproc/gene_expression/Final_data/transcriptome_subcohort_f1_residuals3_v2.Rdata")
eset <- transcriptome_subcohort_f1_residuals3
feat <- featureData(eset)
anno <- pData(feat)
names(anno)
nrow(anno)

# Annotated
an <- anno[anno$GeneSymbol != "",]
nrow(an)
# Non-annotated
nan <- anno[anno$GeneSymbol == "",]
nrow(nan)

# Annotated?
anno$Annotated <- ifelse(anno$GeneSymbol != "", "yes", "no")

# Chr & TCs
chr <- as.data.frame(table(as.character(anno$chromosome)))
write.csv(chr, "chr_ex.csv")

# Locus type
table(as.character(anno$locus.type))

######################################

# Double table (annotated - coding)
table(anno$Annotated, anno$locus.type)

# TC - Amount of genes
tcgenes <- as.character(anno$GeneSymbol)
tcgenes2 <- strsplit(tcgenes, ";")
counts <- lapply(tcgenes2, length)
counts <- unlist(counts)
newdf <- data.frame("TC_id" = anno$transcript_cluster_id, "counts" = counts)
table(counts)
# write.csv(counts, "tc_genecount.csv")

# Genes - N TCs
gene.s <- paste(anno$GeneSymbol, collapse = ";")
gene.v <- strsplit(gene.s, ";")[[1]]
genes <- unique(gene.v)

tccount <- function(x, arg1) {
  sum(grepl(x, arg1))
}

tcc <- lapply(genes, tccount, arg1 = anno$GeneSymbol)
tcc2 <- unlist(tcc)
genes_tcc <- data.frame("Genes" = genes, "TC_count" = tcc2)
# write.csv(genes_tcc, "genes_tcc.csv")

# Genes in > 1 TC
# genes_tcc2 <- genes_tcc[genes_tcc$TC_count > 1,]
# write.csv(genes_tcc2, "genes_tcc2.csv")

# TCs with > 1 Gene
# tc2g <- anno[grep(";",anno$GeneSymbol),]
# gs <- as.character(tc2g$GeneSymbol)
# gs <- strsplit(gs, ";")
# counts <- lapply(gs, length)
# counts <- unlist(counts)
# dftc <- data.frame("TC" = rownames(tc2g), "Genes" = tc2g$GeneSymbol,
#                    "Count" = counts, "Chromosome" = tc2g$chromosome)
# write.csv(dftc, "tc_genes2.csv")
