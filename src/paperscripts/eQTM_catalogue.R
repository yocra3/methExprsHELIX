#'##############################################################################
#' Generate eQTM catalogue figures and tables
#' This script contains the code for
#'  - Volcano plot
#'  - QC checks
#'  - Get number significants pairs and make subgroups
#'  - Make distance plots
#'  - Compare Illumina annotation vs eQTMs
#'  - Compare Non-cell adjusted  (reference) vs cell adjusted models
#'  - Compare with other eQTMs studies
#'##############################################################################

## Load libraries ####
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(S4Vectors)
library(MultiDataSet)
library(ggforce)
library(topGO)
library(FlowSorted.Blood.450k)

## Load datasets ####
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/gexpAnnotation.Rdata")

## Change name when loading results
load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df
featsU <- featStatsDF

load("results/MethComBatExpResidualsCellAdj/allres_simP_cpgs.Rdata")
modC <- df
featsC <- featStatsDF

## Create useful vars
### data.frame only with significant pairs
sigDf <- modU %>%
  as_tibble() %>%
  filter(sigPair)

codingTCs <- subset(expAnnot, Coding == "coding")$transcript_cluster_id
sigTCs <- unique(sigDf$TC)

### data.frame with results and overlaps info
modU_comp <- as_tibble(inner_join(modU, overDF, by = c("CpG", "TC")))

### Modify methylation annotation
methyAnnot <- methyAnnot %>%
  as_tibble() %>%
  mutate(CpG = Name)

# General Overview ####
## Statistics ####
### Significant pairs
nrow(sigDf)
# [1] 63831

### Significant TC
length(unique(sigDf$TC))
# [1] 11071

### Significant coding TC
sum(unique(sigDf$TC) %in% codingTCs)
# [1] 7874

### Significant CpGs
length(unique(sigDf$CpG))
# [1] 35228

## Distribution CpGs/TC
CpG_plot <- sigDf %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = n)) + geom_histogram(binwidth = 1) +
  scale_x_continuous("", limits = c(0, 50)) +
  scale_y_continuous("")  +
  ggtitle("TCs associated with each CpG")

TC_plot <- sigDf %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = n)) + geom_histogram(binwidth = 1) +
  scale_x_continuous("", limits = c(0, 50)) +
  scale_y_continuous("")  +
  ggtitle("CpGs associated with each TC")

png("paper/eQTMs_CpGs_TC_distr.png", width = 2000, height = 1500, res = 300)
plot_grid(CpG_plot, TC_plot, labels = "AUTO", nrow = 2)
dev.off()

sigDf %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   2.000   5.766   6.000 176.000

sigDf %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   1.000   1.812   2.000  26.000

## Distribution pairs per chromosomes ####
png("paper/eQTMs_Chr_distr.png", width = 2000, height = 1500, res = 300)
sigDf %>% 
  mutate(chr = substring(TC, 3, 4),
         chr = gsub("^0", "", chr)) %>%
  group_by(chr) %>%
  summarize(n = n()) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X"))) %>%
  ggplot(aes(x = chr, y = n)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Chromosome") +
  scale_y_continuous(name = "Number Pairs") +
  theme_bw()
dev.off()


## Volcano plot ####
png("paper/volcano_all.png", width = 2000, height = 2000, res = 300)
volcano_plot(modU$p.value, modU$FC/100, paste(modU$CpG, modU$TC), 
             tPV = -log10(1e-8), tFC = 0.01, show.labels = FALSE)
dev.off()


## QCs ####
### CpGs variability
featsU_var <- methyAnnot %>%
  as_tibble() %>%
  mutate(feat = Row.names) %>%
  select(feat, meth_range, variability) %>%
  right_join(featsU, by = "feat") %>%
  mutate(sig = ifelse(p.val.adj < 0.05, "significant", "random"))
table(featsU_var$variability, featsU_var$sig)
#           random significant
# invariant 169446        4647
# variant   181744       30581
chisq.test(table(featsU_var$variability, featsU_var$sig))
# data:  table(featsU_var$variability, featsU_var$sig)
# X-squared = 15894, df = 1, p-value < 2.2e-16

### TC call rate
int_TC <- expAnnot %>%
  select(probeset_id, Expressed, CallRate) %>%
  mutate(sig = ifelse(probeset_id %in% sigTCs, "significant", "random"))
table(int_TC$Expressed, int_TC$sig)
#               random significant
# Expressed      42200       10733
# Not-Expressed   7421         338

chisq.test(table(int_TC$Expressed, int_TC$sig))
# X-squared = 1149, df = 1, p-value < 2.2e-16


### Classify CpG in groups ####
CpGsSum <- modU %>%
  group_by(CpG) %>%
  summarise(Type = ifelse(sum(sigPair) == 0, "Non-significant",
                          ifelse(sum(sigPair) == 1, "Mono", "Multi")),
            Direction = ifelse(sum(sigPair) == 0, "Non-significant",
                               ifelse(all(FC[sigPair] > 0), "Positive", 
                                      ifelse(all(FC[sigPair] < 0), "Inverse", "Both")))) %>%
  mutate(Combined = ifelse(Type == "Non-significant", 
                           "Non-significant", 
                           paste(Type, Direction, sep = "_")))
t <- table(CpGsSum$Type, CpGsSum$Direction)
t <- t[c("Mono", "Multi"), c("Inverse", "Positive", "Both")]
t
#       Inverse Positive  Both
# Mono    12637     8961     0
# Multi    6135     3530  3965
chisq.test(t[, -3])
# X-squared = 68.441, df = 1, p-value < 2.2e-16

sink("paper/CpGs_type.txt")
addmargins(t)
sink()

# Distance plots ####
# Distance distribution
## Signif vs no-signif
dist_all <- ggplot(modU_comp, aes(x = Distance, color = sigPair)) + geom_density() + 
  theme_bw() + 
  scale_color_discrete(name = "", labels = c("Non-significant", "Significant")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "", breaks = NULL) + 
  ggtitle("CpG-TC Distance (all pairs)") + 
  theme(plot.title = element_text(hjust = 0.5))

## Signif per groups
dist_groups <- modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  ggplot(aes(x = Distance, color = Direction)) + geom_density() + 
  theme_bw() + 
  scale_color_manual(name = "", 
                     breaks = c("Inverse", "Positive", "Both"),
                     values = c("grey30", "red", "blue")) +
  facet_grid(~ Type) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "", breaks = NULL) + 
  ggtitle("CpG-TC Distance (Significant pairs)") + 
  theme(plot.title = element_text(hjust = 0.5))

png("paper/distance_distribution.png", width = 2000, height = 2000, res = 300)
plot_grid(dist_all, dist_groups, nrow = 2)
dev.off()

# Distance vs Effect size
png("paper/distance_effect.png", width = 2000, height = 1300, res = 300)
modU_comp %>%
  filter(sigPair) %>%
  mutate(Direction = ifelse(FC > 0, "Positive", "Inverse")) %>%
  ggplot(aes(x = Distance, y = FC/100, color = Direction)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(name = "", 
                     breaks = c("Inverse", "Positive"),
                     values = c("red", "blue")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "FC per 1% Methylation") + 
  ggtitle("CpG-TC Distance vs Effect") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("paper/effect_distribution.png", width = 2000, height = 1500, res = 300)
modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  ggplot(aes(x = abs(FC/100), color = Direction)) + 
  geom_density() +
  theme_bw() + 
  scale_color_manual(name = "", 
                     breaks = c("Inverse", "Positive", "Both"),
                     values = c("grey30", "red", "blue")) +
  facet_grid(~ Type) +
  scale_x_continuous(name = "abs FC per 1% Methylation", 
                     limits = c(0, 0.2)) +
  scale_y_continuous(name = "", breaks = NULL) + 
  ggtitle("Effect size Distribution (Significant pairs)") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Illumina Annotation ####
modU_Annot <- modU %>%
  as_tibble() %>%
  select(CpG, TC, FC, sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  left_join(select(methyAnnot, CpG, UCSC_RefGene_Name), by = "CpG") %>%
  left_join(select(expAnnot, TC, GeneSymbol_Affy), by = "TC") %>%
  mutate(GeneAffy = strsplit(GeneSymbol_Affy, ";"))

methGenes <- unique(unlist(modU_Annot$UCSC_RefGene_Name))
expGenes <- unique(unlist(modU_Annot$GeneAffy))
comGenes <- intersect(methGenes, expGenes)

modU_Annot_f <- modU_Annot %>%
  filter(sigPair)

methGenesf <- unique(unlist(modU_Annot_f$UCSC_RefGene_Name))
expGenesf <- unique(unlist(modU_Annot_f$GeneAffy))
comGenesf <- intersect(methGenesf, expGenesf)

modU_Annot_f %>%
  select(CpG, UCSC_RefGene_Name) %>%
  distinct() %>%
  summarize(n = n(),
            n_g = sum(UCSC_RefGene_Name != ""),
            n_GE = sum(sapply(UCSC_RefGene_Name, function(x) any(x %in% comGenes))))


modU_Annot_f2 <- modU_Annot_f %>%
  filter(sapply(UCSC_RefGene_Name, function(x) any(x %in% comGenes))) %>%
  mutate(match = GeneAffy %in% UCSC_RefGene_Name)


table(modU_Annot_f2$match)
modU_Annot_f2 %>%
  mutate(Dir = ifelse(FC > 0, "Positive", "Inverse")) %>%
  group_by(Dir) %>%
  summarize(N = sum(match),
            P = mean(match))



annotated <- modU_Annot_f2 %>%
  group_by(CpG) %>%
  summarize(found = any(match)) %>%
  left_join(CpGsSum, by = "CpG")
mean(annotated$found)
table(annotated$found)


annotated %>%
  group_by(Combined) %>%
  summarize(N = sum(found),
            P = mean(found))

# Compare models ####
sigDf2 <- modC %>%
  as_tibble() %>%
  filter(sigPair)
## Statistics ####
### Significant pairs
nrow(sigDf2)
# [1] 39749

### Significant TC
length(unique(sigDf2$TC))
# [1] 8886

### Significant coding TC
sum(unique(sigDf2$TC) %in% codingTCs)
# [1] 6288

### Significant CpGs
length(unique(sigDf2$CpG))
# [1] 21966

sigDf2 %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   2.000   4.473   5.000 128.000

sigDf2 %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.00    1.00    1.00    1.81    2.00   24.00


### Classify CpG in groups ####
CpGsSum2 <- modC %>%
  group_by(CpG) %>%
  summarise(Type = ifelse(sum(sigPair) == 0, "Non-significant",
                          ifelse(sum(sigPair) == 1, "Mono", "Multi")),
            Direction = ifelse(sum(sigPair) == 0, "Non-significant",
                               ifelse(all(FC[sigPair] > 0), "Positive", 
                                      ifelse(all(FC[sigPair] < 0), "Inverse", "Both")))) %>%
  mutate(Combined = ifelse(Type == "Non-significant", 
                           "Non-significant", 
                           paste(Type, Direction, sep = "_")))
t <- table(CpGsSum2$Type, CpGsSum2$Direction)
t <- t[c("Mono", "Multi"), c("Inverse", "Positive", "Both")]
t


#       Inverse Positive Both
# Mono     8084     5681    0
# Multi    3738     2400 2063
chisq.test(t[, -3])
# X-squared = 8.2039, df = 1, p-value = 0.00418

sink("paper/CpGs_type2.txt")
addmargins(t)
sink()

## Data.frame with groups comparison
### General
adjList <- list(pairs = paste(sigDf$CpG, sigDf$TC, sep = "-"),
                 TCs = unique(sigDf$TC),
                 CpG = unique(sigDf$CpG))
adjList$codTC <- adjList$TCs[adjList$TCs %in% codingTCs]
cellList <- list(pairs = paste(sigDf2$CpG, sigDf2$TC, sep = "-"),
                 TCs = unique(sigDf2$TC),
                 CpG = unique(sigDf2$CpG))
cellList$codTC <- cellList$TCs[cellList$TCs %in% codingTCs]


modCompDF <- data.frame(Adj = lengths(adjList),
                    Cell = lengths(cellList),
                    Shared = unlist(Map(function(x, y) length(intersect(x, y)), adjList, cellList)),
                    AdjSp = unlist(Map(function(x, y) length(setdiff(x, y)), adjList, cellList)),
                    CellSp = unlist(Map(function(x, y) length(setdiff(x, y)), cellList, adjList)))
sink("paper/ModelsStatisticComparison.txt")
modCompDF[c(1, 2, 4, 3), ]
sink()

### Per CpG type
adjCpgs <- CpGsSum %>% 
  group_by(Combined) %>%
  group_split()
names(adjCpgs) <- sapply(adjCpgs, function(x) x$Combined[1])
adjCpgs <- sapply(adjCpgs, pull, CpG)

cellCpgs <- CpGsSum2 %>% 
  group_by(Combined) %>%
  group_split()
names(cellCpgs) <- sapply(cellCpgs, function(x) x$Combined[1])
cellCpgs <- sapply(cellCpgs, pull, CpG)

modCompDF2 <- data.frame(Adj = lengths(adjCpgs),
                        Cell = lengths(cellCpgs),
                        Shared = unlist(Map(function(x, y) length(intersect(x, y)), adjCpgs, cellCpgs)),
                        AdjSp = unlist(Map(function(x, y) length(setdiff(x, y)), adjCpgs, cellCpgs)),
                        CellSp = unlist(Map(function(x, y) length(setdiff(x, y)), cellCpgs, adjCpgs)))
sink("paper/ModelsStatisticComparisonCpG.txt")
modCompDF2[c(1, 2, 4, 5, 3, 6), ]
sink()

parPlot <- CpGsSum %>%
  left_join(select(CpGsSum2, CpG, Combined), by = "CpG") %>%
  group_by(Combined.x, Combined.y) %>%
  summarize(n = n()) %>%
  filter(!(Combined.x == "Non-significant" & Combined.y == "Non-significant")) %>%
  ungroup() %>%
  mutate(Combined.x = factor(Combined.x, 
                             levels = c("Mono_Inverse", "Mono_Positive", "Multi_Inverse", 
                                        "Multi_Positive", "Multi_Both", "Non-significant")),
         Combined.y = factor(Combined.y, 
                             levels = c("Mono_Inverse", "Mono_Positive", "Multi_Inverse", 
                                        "Multi_Positive", "Multi_Both", "Non-significant"))) %>%
  gather_set_data(1:2) %>%
  ggplot(aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = Combined.y), alpha = 0.6, axis.width = 0.4) +
    geom_parallel_sets_axes(axis.width = 0.4, fill = "gray95",
                            color = "gray80") +
    geom_parallel_sets_labels(colour = 'gray35', size = 4.5, angle = 0, fontface = "bold") +
  scale_x_discrete(labels = c("Adjusted", "Cell Adjusted")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.x  = element_blank()
  )

png("paper/CompModelsCpGs.png", width = 2500, height = 2500, res = 300)
parPlot
dev.off()

## Compare estimates ####
### Merge dataset
mergeTB <- modU %>%
  left_join(modC, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  filter(sigPair.x == TRUE | sigPair.y == TRUE) %>%
  mutate(sigType = ifelse(sigPair.x == TRUE, ifelse(sigPair.y == TRUE, "Both", "Adjusted"), "Cell"))


png("paper/CompModelsP_values.png", width = 2500, height = 2500, res = 300)
ggplot(mergeTB, aes(x = -log10(p.value.y), y = -log10(p.value.x), col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Adjusted") + 
  ggtitle("-log10 p-values comparative") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "")
dev.off()

# mergeTB <- mergeTB %>%
#   mutate(Diff.pval = pnorm(abs(FC.x - FC.y)/sqrt(SD.x**2 + SD.y**2), lower.tail = FALSE), 
#          isDiff = ifelse(Diff.pval < 0.05, "Different", "Equal"))

png("paper/CompModelsEstimates.png", width = 2500, height = 2500, res = 300)
ggplot(mergeTB, aes(x = FC.y, y = FC.x, col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Adjusted") + 
  ggtitle("Estimates comparative") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "")
dev.off()


## Run enrichment on specific TCs ####
### Define function to run enrichment
computeGOs <- function(df){
  
  genesv <- df$sig
  names(genesv) <- df$TC
  
  Data <- new("topGOdata", 
              description = "GO analysis of TCs",
              ontology = "BP",
              allGenes = genesv,
              annot = annFUN.db,
              nodeSize = 10,
              affyLib = "hta20transcriptcluster.db")
  
  class <- runTest(Data, algorithm = "classic", statistic = "fisher")
  elim <- runTest(Data, algorithm = "elim", statistic = "fisher")
  weigth <- runTest(Data, algorithm = "weight", statistic = "fisher")
  mixed <- runTest(Data, algorithm = "weight01", statistic = "fisher")
  
  
  finTab <- GenTable(Data, 
                     classic = class, 
                     elim = elim, 
                     weight = weigth, 
                     w0 = mixed,
                     orderBy = "weight",
                     ranksOf = "elim",
                     topNodes = length(score(class)))
  gos <- list(classic = class, elim = elim,
              weigth = weigth, 
              mixed = mixed)
  list(gos = gos, table = finTab)
}
adjTCs <- setdiff(adjList$TC, cellList$TC)
adjGenes <- modU %>%
  group_by(TC) %>%
  summarize(sig = factor(
    ifelse(any(TC %in% adjTCs & sigPair), 1, 0)))
adjGOs <- computeGOs(adjGenes)

a <- adjGOs$table %>%
  filter(w0 < 0.01 & classic < 0.01) %>%
  arrange(w0)
write.table(a[, c("GO.ID", "Term", "w0", "classic")], file = "paper/GOsAdjSp.txt", 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)

cellTCs <- setdiff(cellList$TC, adjList$TC)
cellGenes <- modC %>%
  group_by(TC) %>%
  summarize(sig = factor(
    ifelse(any(TC %in% cellTCs & sigPair), 1, 0)))
cellGOs <- computeGOs(cellGenes)

a <- cellGOs$table %>%
  filter(w0 < 0.01 & classic < 0.01) %>%
  arrange(w0)
write.table(a[, c("GO.ID", "Term", "w0", "classic")], file = "paper/GOsCellSp.txt", 
            quote = FALSE, col.names = TRUE,
            sep = "\t", row.names = FALSE)


## Compare CpGs with Reinius ####
data(FlowSorted.Blood.450k.compTable)
adjCpGs <- setdiff(adjList$CpG, cellList$CpG)
cellCpGs <- setdiff(cellList$CpG, adjList$CpG)

data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(adjCpGs, cellCpGs), "Fstat"],
           Type = rep(c("Adj", "Cell"), c(length(adjCpGs), length(cellCpGs)))) %>%
  lm(formula = log10(Fstat) ~ Type, data = .) %>%
  summary()
  
# Call:
#   lm(formula = log10(Fstat) ~ Type, data = .)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -3.6491 -0.3193 -0.0387  0.3510  2.1680
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.26598    0.00460  275.21   <2e-16 ***
#   TypeCell     0.16409    0.01114   14.73   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5934 on 20061 degrees of freedom
# (73 observations deleted due to missingness)
# Multiple R-squared:  0.0107,    Adjusted R-squared:  0.01065
# F-statistic: 216.9 on 1 and 20061 DF,  p-value: < 2.2e-16
# 

png("paper/CompModelsCpGCellSpecific.png", width = 1500, height = 1000, res = 300)
data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(adjCpGs, cellCpGs), "Fstat"],
           Type = rep(c("Adjusted", "Cell"), c(length(adjCpGs), length(cellCpGs)))) %>%
  ggplot(aes(y = log10(Fstat), x = Type, fill = Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "")
dev.off()


## Compare with other eQTM studies ####
### Bonder
### Gutierrez-Arcelus
