###############################################################
# ANNOTATION DESCRIPTIVES 2 (EXPRESSION)
###############################################################

library('Biobase')

load("summ_exp_sm_1.RData")
anno <- rowData(se)
gr <- rowRanges(se)

table(anno$geneNA)

# Chr & TCs
chr <- as.data.frame(table(as.character(seqnames(gr))))
write.csv(chr, "chr_expr2.csv")

######################################

# Double table (annotated - coding)
table(anno$geneNA, anno$phase)

# TC - Amount of genes
tcg.f <- function(x) {
  if (x == "" | is.na(x)) {0}
  else {length(x)}
}
counts.affy <- lapply(anno$GeneSymbol_Affy, tcg.f)
counts.affy <- unlist(counts.affy)
anno$counts.affy <- counts.affy
counts.db <- lapply(anno$GeneSymbolDB2,tcg.f)
counts.db <- unlist(counts.db)
anno$counts.db <- counts.db

# newdf <- anno[, c(1,3,7,10,11,13,14,15,16,17,18)]
newdf <- anno[,c(1,3,7,16,17,18)]
# write.csv(newdf, "tc_genecount2.csv")
as.data.frame(table(counts.affy))
as.data.frame(table(counts.db))

# Genes - N TCs
gs1 <- unlist(anno$GeneSymbol_Affy)
gs1 <- unique(gs1)
gs2 <- unlist(anno$GeneSymbolDB2)
gs2 <- unique(gs2)

tccount <- function(x, arg1) {
  x <- paste("^", x, "$", sep = "")
  sum(grepl(x, arg1))
}

tcc1 <- lapply(gs1, tccount, arg1 = anno$GeneSymbol_Affy)
tcc1 <- unlist(tcc1)
tcc2 <- lapply(gs2, tccount, arg1 = anno$GeneSymbolDB2)
tcc2 <- unlist(tcc2)

genes_tcc1 <- data.frame("Genes" = gs1, "TC_count" = tcc1)
genes_tcc1 <- genes_tcc1[order(genes_tcc1$TC_count),]
genes_tcc2 <- data.frame("Genes" = gs2, "TC_count" = tcc2)
genes_tcc2 <- genes_tcc2[order(genes_tcc2$TC_count),] 
# write.csv(genes_tcc1, "genes_tcc1.csv")
# write.csv(genes_tcc2, "genes_tcc2.csv")