###############################################################################
# Reestimate pvalues using simulations
###############################################################################
## Run in ~/data/smari/
library(parallel)

## Load original results
load("data_res2/all_output.RData")

## Load list of folders
folders <- dir("simulations", pattern = "sim[0-9]")
chr <- 1:22
  
pvalsGenes <- lapply(chr, function(chrom){
  print(chrom)
  resList <- lapply(folders, function(simFold) 
    read.table(paste0("simulations/", simFold, "/output", chrom, ".txt"), as.is = TRUE))
  
  genes <- unique(resList[[1]]$V2)
  
  res <- mclapply(genes, function(gene){
    resGene <- lapply(resList, function(x) {
      df <- x[x$V2 == gene, ]
      pvals <- df$V5
      })
    perms <- Reduce(cbind, resGene)
    pmin <- apply(perms, 2, min, na.rm=TRUE)
    
    llhd2 <- function(x,p) {
      ans <- -sum(log(dbeta(x,p[1],p[2])))
      ans
    }
    
    pIni <- c(1,100)
    
    param <- nlm(llhd2, x=pmin, p=pIni)$estimate
    
    minObs <- min(all_output[all_output$TC == gene, "p.value"])
    pbeta(minObs, param[1], param[2], lower.tail = TRUE)
  }, mc.cores = 5)
  res <- unlist(res)
  names(res) <- genes
  res
})
newPvals <- unlist(pvalsGenes)
selGenes <- names(newPvals[p.adjust(newPvals, method = "bonferroni") < 0.05])

pvalsOld <- tapply(all_output$p.value, all_output$TC, min)

## Select threshold
max(pvalsOld[selGenes])
# [1] 9.024952e-07

selPairs <- subset(all_output, p.value <= max(pvalsOld[selGenes]) & TC %in% selGenes)

## Add annotation
load("overlaps_d.RData")

rownames(selPairs) <- paste(selPairs$CpG, selPairs$TC, sep = "-")
rownames(overlaps_d) <- paste(overlaps_d$CpG, overlaps_d$TC, sep = "-")

selPairs$dist <- overlaps_d[rownames(selPairs), "dif"]

save(selPairs, file = "simulations/selectedPairs.Rdata")

