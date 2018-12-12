#'##############################################################################
#' Compare models adjusted and not adjusted
#'##############################################################################

# Load libraries ####
library(MEAL)
library(minfi)
library(limma)

# Load data ####
load("data_res2/all_output.RData")
mod2 <- all_output

load("data_res4/all_output.Rdata")
mod4 <- all_output

## Methylation data
load("gset_sm_1.RData")

## Expression data
load("summ_exp_sm_res2.RData")

## Associate CpGs to cell counts ####
model <- model.matrix(~NK + Bcell + CD4T + CD8T + Eos + Neu + 
                        Mono + cohort + e3_sex + age_sample_years,
                      data = colData(gset))
lm <- lmFit(getBeta(gset), model)
lm <- eBayes(lm)

## Associate Genes to cell counts ####
modelExp <- model.matrix(~NK + Bcell + CD4T + CD8T + Eos + Neu + 
                        Mono + cohort + e3_sex + age_sample_years,
                      data = colData(se_res))
lmExp <- lmFit(assay(se_res), model)
lmExp <- eBayes(lmExp)



## model A: not including cell counts
## model C: including cell counts

## Select pairs exclusive of model A ####
crudepairs <- -log10(mod2$p.value) < 8 & -log10(mod4$p.value) > 8
mod4sel <- mod4[crudepairs, ]

## Make summary data.frame
crudeDF <- mod4sel[, 1:5]
crudeDF <- cbind(crudeDF, mod2[crudepairs, 3:5])
crudeDF <- cbind(crudeDF, data.frame(rowRanges(gset)[crudeDF$CpG])[, c("UCSC_RefGene_Name", "UCSC_RefGene_Group", "chr", "pos")])
crudeDF <- cbind(crudeDF, rowData(se_res[crudeDF$TC,])[, c("TSS_Affy", "GeneSymbol_Affy")])

crudeDF <- crudeDF[, c(1, 9:12, 2, 14, 13, 3:8)]
colnames(crudeDF) <- c("CpG", "CpG-Gene", "CpG-Gene Position", "CpG-Chr", "CpG-Pos", 
                       "TC", "TC-Gene", "TC-TSS", "eQTM-Estimate_Crude", "eQTM-SE_Crude", 
                       "eQTM-p.value_Crude", "eQTM-Estimate_Cells", "eQTM-SE_Cells",
                       "eQTM-p.value_Cells")

crudeDF[["CpG-Gene"]] <- sapply(crudeDF[["CpG-Gene"]], paste, collapse = ";")
crudeDF[["TC-Gene"]] <- sapply(crudeDF[["TC-Gene"]], paste, collapse = ";")

crudeDF[["eQTM-Estimates-Diff_pvalue"]] <- pnorm(abs(
  crudeDF[["eQTM-Estimate_Crude"]] - crudeDF[["eQTM-Estimate_Cells"]])/
    sqrt(crudeDF[["eQTM-SE_Crude"]]**2 + crudeDF[["eQTM-SE_Cells"]]**2), 
  lower.tail = FALSE)

## Add cell counts association
cpgCells <- lapply(2:8, function(coef) {
  df <- topTable(lm, coef = coef, n = Inf)
  df[crudeDF$CpG, 4]
})

ExpCells <- lapply(2:8, function(coef) {
  df <- topTable(lmExp, coef = coef, n = Inf)
  df[crudeDF$TC, 4]
})

crudeDF <- cbind(crudeDF, Reduce(cbind, cpgCells), Reduce(cbind, ExpCells))
colnames(crudeDF)[16:29] <- paste(c("NK", "Bcell", "CD4T", "CD8T", "Eos", "Neu", "Mono"), 
                                  rep(c("CpG", "TC"), each = 7), 
                                  "p.value", sep = "-")
write.table(crudeDF, file = "nonCellAdjustedResults.txt", col.names = TRUE, 
            row.names = FALSE, quote = FALSE)

## Select pairs exclusive of model C ####
cellpairs <- -log10(mod2$p.value) > 8 & -log10(mod4$p.value) < 8
mod2sel <- mod2[cellpairs, ]

## Make summary data.frame
cellDF <- mod2sel[, 1:5]
cellDF <- cbind(cellDF, mod4[cellpairs, 3:5])
cellDF <- cbind(cellDF, data.frame(rowRanges(gset)[cellDF$CpG])[, c("UCSC_RefGene_Name", "UCSC_RefGene_Group", "chr", "pos")])
cellDF <- cbind(cellDF, rowData(se_res[cellDF$TC,])[, c("TSS_Affy", "GeneSymbol_Affy")])

cellDF <- cellDF[, c(1, 9:12, 2, 14, 13, 3:8)]
colnames(cellDF) <- c("CpG", "CpG-Gene", "CpG-Gene Position", "CpG-Chr", "CpG-Pos", 
                       "TC", "TC-Gene", "TC-TSS", "eQTM-Estimate_Cells", "eQTM-SE_Cells", 
                       "eQTM-p.value_Cells", "eQTM-Estimate_Crude", "eQTM-SE_Crude",
                       "eQTM-p.value_Crude")

cellDF[["CpG-Gene"]] <- sapply(cellDF[["CpG-Gene"]], paste, collapse = ";")
cellDF[["TC-Gene"]] <- sapply(cellDF[["TC-Gene"]], paste, collapse = ";")


cellDF[["eQTM-Estimates-Diff_pvalue"]] <- pnorm(abs(
  cellDF[["eQTM-Estimate_Cells"]] - cellDF[["eQTM-Estimate_Crude"]])/
    sqrt(cellDF[["eQTM-SE_Crude"]]**2 + cellDF[["eQTM-SE_Cells"]]**2), 
  lower.tail = FALSE)


## Add cell counts association
cpgCells <- lapply(2:8, function(coef) {
  df <- topTable(lm, coef = coef, n = Inf)
  df[cellDF$CpG, 4]
})

ExpCells <- lapply(2:8, function(coef) {
  df <- topTable(lmExp, coef = coef, n = Inf)
  df[cellDF$TC, 4]
})


cellDF <- cbind(cellDF, Reduce(cbind, cpgCells), Reduce(cbind, ExpCells))
colnames(cellDF)[16:29] <- paste(c("NK", "Bcell", "CD4T", "CD8T", "Eos", "Neu", "Mono"), 
                                  rep(c("CpG", "TC"), each = 7), 
                                  "p.value", sep = "-")
write.table(cellDF, file = "cellAdjustedResults.txt", col.names = TRUE, 
            row.names = FALSE, quote = FALSE)
