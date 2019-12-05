#'##############################################################################
#' Generate eQTM interpretation figures and tables
#' This script contains the code for
#'  - Gene enrichment
#'  - CpG enrichment 
#'  - eQTMs and genes and CpGs changing during childhood
#'  - Overlap with EWAS atlas
#'##############################################################################

# Load data and libraries ####
library(topGO)
library(ggplot2)
library(cowplot)
library(tidyr)
library(hta20transcriptcluster.db)
library(FlowSorted.Blood.450k)
library(openxlsx)
library(GenomicRanges)
library(dplyr)

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/gexpAnnotation.Rdata")

# Get useful variables ####
CpGsSum <- df %>%
  group_by(CpG) %>%
  summarise(Type = ifelse(sum(sigPair) == 0, "Non-significant",
                          ifelse(sum(sigPair) == 1, "Mono", "Multi")),
            Direction = ifelse(sum(sigPair) == 0, "Non-significant",
                               ifelse(all(FC[sigPair] > 0), "Positive", 
                                      ifelse(all(FC[sigPair] < 0), "Inverse", "Both")))) %>%
  mutate(Combined = ifelse(Type == "Non-significant", 
                           "Non-significant", 
                           paste(Type, Direction, sep = "_")))



## Gene enrichment ####
### Define functions
#' Compute GOs compute GO using two tests: 
#'  - classic: get p-value for each GO independent of the other
#'  - weight01: takes into account the GO hierarchy
## We selected only those terms with a p-value < 0.01 in both methods.

getGOdata <- function(df){
  genesv <- df$sig
  names(genesv) <- df$TC
  
  Data <- new("topGOdata", 
              description = "GO analysis of TCs",
              ontology = "BP",
              allGenes = genesv,
              annot = annFUN.db,
              nodeSize = 10,
              affyLib = "hta20transcriptcluster.db")
}

computeGOs <- function(df){
  
  Data <- getGOdata(df)
  mixed <- runTest(Data, algorithm = "weight01", statistic = "fisher")
  
  
  finTab <- GenTable(Data, 
                     w0 = mixed,
                     orderBy = "w0",
                     topNodes = length(score(mixed)))
  list(go = mixed, table = finTab)
}

runGO <- function(annot, df, filter){
  selCpGs <- filter(annot, Combined == filter)$CpG
  selGenes <- df %>%
    group_by(TC) %>%
    summarize(sig = factor(
      ifelse(any(CpG %in% selCpGs & sigPair), 1, 0)))
  go <- computeGOs(selGenes)
  go
}


## All genes
allGenes <- df %>%
  group_by(TC) %>%
  summarize(sig = factor(ifelse(any(sigPair), 1, 0)))
allGos <- computeGOs(allGenes)
allGOData <- getGOdata(allGenes)

## Subtypes
subtypes <- c("Mono_Inverse", "Mono_Positive", "Multi_Both", "Multi_Inverse", 
              "Multi_Positive")
names(subtypes) <- subtypes
subtypes_go <- lapply(subtypes, runGO, annot = CpGsSum, df = df)

## Save as Rdata data and process in local
save(allGos, subtypes_go, file = "paper/GOobjects.Rdata")

## Genes used in GOs
allGos$go@geneData["Significant"]
# Significant
# 6675

sapply(subtypes_go, function(x) x$go@geneData["Significant"])
# Mono_Inverse.Significant  Mono_Positive.Significant
# 3389                       2637
# Multi_Both.Significant  Multi_Inverse.Significant
# 2591                       2810
# Multi_Positive.Significant
# 1872

## GOs changed
sapply(subtypes_go, function(x) sum(x$table$w0 < 0.001))
# Mono_Inverse  Mono_Positive     Multi_Both  Multi_Inverse Multi_Positive
# 24             18             51             46             42



## It does not work on server. Run locally.
library(GOfuncR)
library(topGO)
library(dplyr)
server <- "//isg10174/data/WS_HELIX/HELIX_analyses/expr_met_SM/paper/"

load(paste0(server, "GOobjects.Rdata"))

## Define categories for GO terms
## Use http://amigo.geneontology.org/amigo/dd_browse and only select GOs > 10 genes in H. Sapiens.
## Exclude terms with regulation
##' Immune: 
##' - immune system process
##' - cell killing
##' Use level one/two parent
##' 
##' Adaptive immunity:
##' - adaptive immune response
##' - B cell selection
##' - T cell selection
##' - antigen processing and presentation
##' - B killer cell activation involved in immune response
##' - T cell activation involved in immune response
##' - complement activation
##' - complement-dependent cytotoxicity
##' - T cell mediated cytotoxicity
##' - B cell mediated immunity
##' - T cell mediated immunity
##' - opsonization
##' - B cell activation involved in immune response
##' - T cell activation involved in immune response
##' - B cell activation
##' - T cell activation
##' - B cell activation involved in immune response
##' - T cell activation involved in immune response
##' - B cell differentiation
##' - T cell differentiation
##' - B cell proliferation
##' - T cell proliferation
##' - lymphocyte homeostasis
##' - T cell extravasation
##' - lymphocyte migration
##' - lymphocyte coestimulation
##' - immunoglobulin production
##' - T cell cytokine production
##' - humoral immune response mediated by circulating immunoglobulin
##' - antigen receptor-mediated signaling pathway
##' 
##' 
##' Innate immunity:
##' - innate immune response
##' - astrocyte activation involved in immune response
##' - endothelial cell activation involved in immune response
##' - myeloid cell activation involved in immune response
##' - natural killer cell activation involved in immune response
##' - natural killer cell mediated cytotoxicity
##' - neutrophil mediated cytotoxicity
##' - immune response-regulating cell surface receptor signaling pathway involved in phagocytosis
##' - leukocyte degranulation
##' - dendritic cell cytokine production
##' - myeloid leukocyte mediated immunity
##' - natural killer cell mediated immunity
##' - myeloid cell activation involved in immune response
##' - natural killer cell activation involved in immune response
##' - leukocyte activation involved in inflammatory response
##' - natural killer cell activation involved in immune response
##' - natural killer cell differentiation
##' - natural killer cell proliferation
##' - natural killer cell activation
##' - myeloid leukocyte activation
##' - neutrophil homeostasis
##' - nuetrophil extravasation
##' - dendritic cell migration
##' - mononuclear cell migration
##' - myeloid leukocyte migration
##' - myeloid cell homeostasis
##' - dendritic cell cytokine production
##' - myeloid leukocyte cytokine production
##' - antimicrobial humoral response
##' - inflammatory response to antigenic stimulus

immune <- c("immune system process", "cell killing")
adaptive <- c("adaptive immune response",
              "antigen receptor-mediated signaling pathway",
              "B cell selection",
              "T cell selection",
              "antigen processing and presentation",
              "B killer cell activation involved in immune response",
              "T cell activation involved in immune response",
              "complement activation",
              "complement-dependent cytotoxicity",
              "T cell mediated cytotoxicity",
              "B cell mediated immunity",
              "T cell mediated immunity",
              "opsonization",
              "B cell activation involved in immune response",
              "T cell activation involved in immune response",
              "B cell activation",
              "T cell activation",
              "B cell activation involved in immune response",
              "T cell activation involved in immune response",
              "B cell differentiation",
              "T cell differentiation",
              "B cell proliferation",
              "T cell proliferation",
              "lymphocyte homeostasis",
              "T cell extravasation",
              "lymphocyte migration",
              "lymphocyte coestimulation",
              "immunoglobulin production",
              "T cell cytokine production",
              "humoral immune response mediated by circulating immunoglobulin")
innate <- c("innate immune response",
            "astrocyte activation involved in immune response",
            "endothelial cell activation involved in immune response",
            "myeloid cell activation involved in immune response",
            "natural killer cell activation involved in immune response",
            "natural killer cell mediated cytotoxicity",
            "neutrophil mediated cytotoxicity",
            "immune response-regulating cell surface receptor signaling pathway involved in phagocytosis",
            "leukocyte degranulation",
            "dendritic cell cytokine production",
            "myeloid leukocyte mediated immunity",
            "natural killer cell mediated immunity",
            "myeloid cell activation involved in immune response",
            "natural killer cell activation involved in immune response",
            "leukocyte activation involved in inflammatory response",
            "natural killer cell activation involved in immune response",
            "natural killer cell differentiation",
            "natural killer cell proliferation",
            "natural killer cell activation",
            "myeloid leukocyte activation",
            "neutrophil homeostasis",
            "nuetrophil extravasation",
            "dendritic cell migration",
            "mononuclear cell migration",
            "myeloid leukocyte migration",
            "myeloid cell homeostasis",
            "dendritic cell cytokine production",
            "myeloid leukocyte cytokine production",
            "antimicrobial humoral response",
            "inflammatory response to antigenic stimulus")

top_GOs <- subset(get_child_nodes("GO:0008150"), distance == 1)$child_name
celProc_GOs <- subset(get_child_nodes("GO:0009987"), distance == 1)$child_name

sumF <- function(x, topGOs){
  selGOs <- topGOs[topGOs %in% x]
  paste(selGOs, collapse = ";")
}


addImmunityInfo <- function(tab){
  get_parent_nodes(tab$GO.ID) %>%
    group_by(child_go_id) %>%
    summarize(parent = ifelse(any(immune %in% parent_name), 
                              "immune", 
                              sumF(parent_name, topGOs = top_GOs)), 
              immune = ifelse(any(adaptive %in% parent_name), 
                              "adaptive",
                              ifelse(any(innate %in% parent_name), 
                                     "innate",
                                     "undetermined")),
              GO_term = head(parent_name, 1)) %>%
    mutate(immune = ifelse(immune == "undetermined" & parent != "immune",
                           NA, immune),
           GO.ID = child_go_id) %>%
    right_join(tab)  %>%
    select(GO.ID, GO_term, w0, parent, immune)
}

## Select GOs with p-value < 0.001
allMod <- addImmunityInfo(subset(allGos$table, as.numeric(w0) < 0.001))
subtypesTabs <- lapply(subtypes_go, function(x){
  subset(addImmunityInfo(x$tab), as.numeric(w0) < 0.001)
})
  
  
write.table(allMod[, c("GO.ID", "GO_term", "w0", "parent", "immune")], 
            file = paste0(server, "/GOsAllGenes.txt"), 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)

lapply(names(subtypesTabs), function(x){
  write.table(subtypesTabs[[x]][, c("GO.ID", "GO_term", "w0", "parent", "immune")], 
              file = paste0(server, "/GOs", x, ".txt"), 
              quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
  
})

## Summarize GOs results
### Including all eQTMs
list(c(n = nrow(allMod), imm = sum(allMod$parent == "immune"), 
       immP = round(mean(allMod$parent == "immune")*100, 1)), 
     table(allMod$immune),
     round(prop.table(table(allMod$immune))*100, 1)
)

tail(sort(sapply(top_GOs, function(x) sum(grepl(x, allMod$parent)))))
# response to stimulus 
# 10 
# positive regulation of biological process 
# 11 
# metabolic process 
# 14 
# regulation of biological process 
# 20 
# biological regulation 
# 21 
# cellular process 
# 28 

## Summarize GOs results
### Per CpG type
lapply(subtypesTabs, function(x){
  list(c(n = nrow(x), imm = sum(x$parent == "immune"), 
    immP = round(mean(x$parent == "immune")*100, 1)), 
    table(x$immune),
    round(prop.table(table(x$immune))*100, 1)
  )
})

lapply(subtypesTabs, function(y){
  tail(sort(sapply(top_GOs, function(x) 
    sum(grepl(x, y$parent)))))
})
### Proportion immune vs others
goTab <- sapply(subtypesTabs, function(x) table(x$parent == "immune"))
chisq.test(goTab)

### Proportion adaptive/innate
immuneTab <- sapply(subtypesTabs, function(x) table(x$immune)[c("adaptive", "innate", "undetermined")])
immuneTab[is.na(immuneTab)] <- 0
chisq.test(immuneTab)

  data.frame() %>%
  t() %>%
  data.frame() %>%
  mutate(group = rownames(.),
         group = ifelse(grepl("Mono", group), "Mono", "Multi")) %>%
  group_by(group) %>%
  summarize(adap = sum(adap),
            innate = sum(innate)) 
  chisq.test()

# CpG Enrichment ####
rownames(methyAnnot) <- methyAnnot$Name
methyAnnot$GeneRel <- ifelse(methyAnnot$UCSC_RefGene_Name == "", "Intergenic", "Genic")
types <- c("Both", "Positive", "Inverse")

## Define vars and function
getOR <- function(cols, df){
  t <- data.matrix(df[, cols])
  or <- t[1, 1]/(t[1, 2])/(t[2, 1])*t[2, 2]
  p.val <- chisq.test(t)$p.value
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
}

g <- function(x, gpos, cols = c("in", "ou")){
  sapply(gpos, function(y) getOR(paste0(y, cols), df = x))
}


## Methylation levels ####
methMethLevels <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median_cat) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  group_by(Direction, median_cat) %>%
  summarize(n = n()) %>%
  spread(median_cat, n) %>%
  mutate(Type0 = ifelse(Direction == "Non-significant", "Non-significant", "Significant"), 
         tot = sum(low, medium, high))

cats <- c("low", "medium", "high")

## Define vars and function
getOR2 <- function(col, df){
  t <- data.matrix(df[, c(col, colnames(df)[ncol(df)])])
  t[, 2] <- t[, 2] - t[, 1]
  or <- t[1, 1]/(t[1, 2])/(t[2, 1])*t[2, 2]
  p.val <- chisq.test(t)$p.value
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
}

g2 <- function(x, gpos){
  sapply(gpos, function(y) getOR2(y, df = x))
}

methMethLevels %>%
  mutate(lowp = low/tot,
         medp = medium/tot,
         highp = high/tot) %>%
  gather(Levels, proportion, 7:9)


methMethLevels %>% 
  select(Direction, low, medium, high) %>%
  gather(Levels, vals, 2:4) %>%
  spread(Direction, vals) %>%
  mutate(Significant = Both + Positive + Inverse) %>%
  gather(CpGtype, vals, 2:5) %>%
  mutate(CpGtype = factor(CpGtype, levels = c("Significant", "Non-significant", types))) %>%
  group_by(CpGtype) %>%
  mutate(prop = vals/sum(vals)) %>%
  ggplot(aes(x = Levels, y = prop, fill = CpGtype)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(name = "Methylation Levels", limits = cats) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#000000", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 


allMethLevs <- methMethLevels %>%
  group_by(Type0) %>%
  select(-Direction) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(cats) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:3) %>%
  mutate(Type = "Significant")

typesMethLevs <- lapply(types, function(t){
  methMethLevels %>%
    filter(Direction %in% c(t, "Non-significant")) %>%
    group_by(Direction) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    g2(cats) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:3) %>%
    mutate(Type = t)
})
combMethLevs <- Reduce(rbind, typesMethLevs)
combMethLevs <- rbind(allMethLevs, combMethLevs)

png("paper/CpGEnrichMethLevels.png", width = 3000, height = 2000, res = 300)
combMethLevs %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Inverse", "Positive", "Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Methylation Levels", limits = cats) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 
dev.off()

## Multi vs mono in medium
methMethLevels %>%
  filter(Combined != "Non-significant") %>%
  mutate(group = ifelse(grepl("Mono", Combined), "Mono", "Multi")) %>%
  group_by(group) %>%
  summarize_if(is.numeric, sum) %>%
  getOR2(col = "medium", df = .)
  

## Negative vs positive in low
methMethLevels %>%
  filter(Combined != "Non-significant") %>%
  mutate(group = ifelse(grepl("Positive", Combined), "Positive", "Other")) %>%
  group_by(group) %>%
  summarize_if(is.numeric, sum) %>%
  getOR2(col = "low", df = .)

## Multi positive vs other in high
methMethLevels %>%
  filter(Combined != "Non-significant") %>%
  mutate(group = ifelse(Combined == "Multi_Positive", "Positive", "Other")) %>%
  group_by(group) %>%
  summarize_if(is.numeric, sum) %>%
  getOR2(col = "high", df = .)

## CpG Islands ####
methIsland <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, Relation_to_Island) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  group_by(Direction, Relation_to_Island) %>%
  summarize(n = n()) %>%
  spread(Relation_to_Island, n) %>%
  mutate(Type0 = ifelse(Direction == "Non-significant", "Non-significant", "Significant"), 
      tot = sum(Island, N_Shelf, N_Shore, OpenSea, S_Shelf, S_Shore))

islandStates <- c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "OpenSea")

allIsland <- methIsland %>%
  group_by(Type0) %>%
  select(-Direction) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(islandStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:length(islandStates)) %>%
  mutate(Type = "Significant")

typesIsland <- lapply(types, function(t){
  methIsland %>%
    filter(Direction %in% c(t, "Non-significant")) %>%
    group_by(Direction) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    g2(islandStates) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:length(islandStates)) %>%
    mutate(Type = t)
})
combIsland <- Reduce(rbind, typesIsland)
combIsland <- rbind(allIsland, combIsland)

png("paper/CpGEnrichIsland.png", width = 3000, height = 2000, res = 300)
combIsland %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Inverse", "Positive", "Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Methylation Levels", limits = islandStates) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 
dev.off()

## CpG Island vs methylation levels
png("paper/medianMethvsIsland.png", width = 3000, height = 2000, res = 300)
top <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, Relation_to_Island, median) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Significant = ifelse(Type == "Non-significant", "Non-eQTM", "eQTM")) %>%
  ggplot(aes(x = Relation_to_Island, y = median, fill = Significant)) +
  geom_boxplot() +
  scale_x_discrete(name = "CpG island", limits = islandStates) +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_discrete(name = "") +
  theme_bw()

down <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, Relation_to_Island, median) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Combined = factor(Combined, levels = c("Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both", "Non-significant"))) %>%
  ggplot(aes(x = Relation_to_Island, y = median, fill = Combined)) +
  geom_boxplot() +
  scale_x_discrete(name = "CpG island", limits = islandStates) +
  scale_fill_manual(name = "CpG Type", values = c("#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73", "#555555")) +
  scale_y_continuous(name = "Median methylation") +
  theme_bw()

plot_grid(top, down, nrow = 2)
dev.off()



## Chromatin states ####
chromStates <- c("TssA", "TssAFlnk", "TxFlnk", "TxWk", "Tx", "EnhG", "Enh",
                 "ZNF.Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
                 "ReprPCWk", "Quies")
sum2 <- function(x) sum(!x)

methChromSt <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, eval(chromStates)) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  group_by(Direction) %>%
  summarize_at(chromStates, list(sum = sum, sum2 = sum2)) %>%
  mutate(Type0 = ifelse(Direction == "Non-significant", "Non-significant", "Significant"))


allChromSt <- methChromSt %>%
  group_by(Type0) %>%
  select(ends_with("sum"), ends_with("sum2")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g(chromStates, cols = c("_sum", "_sum2")) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:15) %>%
  mutate(Type = "Significant")

typesChromSt <- lapply(types, function(t){
  methChromSt %>%
    filter(Direction %in% c(t, "Non-significant")) %>%
    group_by(Direction) %>%
    select(ends_with("sum"), ends_with("sum2")) %>%
    summarize_all(list(sum)) %>%
    g(chromStates, cols = c("_sum", "_sum2")) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:15) %>%
    mutate(Type = t)
})
combChromSt <- Reduce(rbind, typesChromSt)
combChromSt <- rbind(allChromSt, combChromSt)

png("paper/CpGEnrichChromStates.png", width = 3000, height = 2000, res = 300)
combChromSt %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, 
                       levels = c("Significant", "Inverse", "Positive", "Both")),
         Group = factor(ifelse(Region %in% c("TssA", "TssAFlnk"), "TssProxProm",
                               ifelse(Region %in% c("Tx", "TxWk"), "ActTrans", 
                                      ifelse(Region %in% c("Enh", "EnhG"), "Enhancer", 
                                             ifelse(Region %in% c("TssBiv", "BivFlnk", "EnhBiv"), "BivReg", 
                                                    ifelse(Region %in% c("ReprPC", "ReprPCWk"), "ReprPoly", Region)
                                                    )
                                             )
                                      )
                               ), 
                        levels = c("TssProxProm", "TxFlnk", "ActTrans", "Enhancer", "ZNF.Rpts", "BivFlnk", "BivReg", "Het", "ReprPoly", "Quies")
                        )
         ) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  facet_wrap(~ Group, scales = "free_x") +
  theme_bw() 
dev.off()


## Chromatin states vs methylation levels
png("paper/medianMethvsChromState.png", width = 3000, height = 2000, res = 300)
  as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, eval(chromStates), median) %>%
  right_join(CpGsSum, by = "CpG") %>%
  gather(Region, Val, 2:16) %>%
  filter(Val == TRUE) %>%
  mutate(Significant = ifelse(Type == "Non-significant", "Non-eQTM", "eQTM"),
         Group = factor(ifelse(Region %in% c("TssA", "TssAFlnk"), "TssProxProm",
                               ifelse(Region %in% c("Tx", "TxWk"), "ActTrans", 
                                      ifelse(Region %in% c("Enh", "EnhG"), "Enhancer", 
                                             ifelse(Region %in% c("TssBiv", "BivFlnk", "EnhBiv"), "BivReg", 
                                                    ifelse(Region %in% c("ReprPC", "ReprPCWk"), "ReprPoly", Region)
                                             )
                                      )
                               )
         ), 
         levels = c("TssProxProm", "TxFlnk", "ActTrans", "Enhancer", "ZNF.Rpts", "BivFlnk", "BivReg", "Het", "ReprPoly", "Quies")
         )
  ) %>%
  ggplot(aes(x = Region, y = median, fill = Significant)) +
  geom_boxplot() +
  facet_wrap(~ Group, scales = "free_x") +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_discrete(name = "") +
  theme_bw()
dev.off()



### Test variables modifying probability of being an eQTM ####
as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Significant = ifelse(Type == "Non-significant", "Non-eQTM", "eQTM"),
         Significant = factor(Significant, levels = c("Non-eQTM", "eQTM")),
         median_trans = -(median-0.5)^2) %>%
  glm(formula(paste("Significant ~ median_trans + Relation_to_Island + GeneRel +",
                    paste(chromStates, collapse = "+"))), 
      family = "binomial", .) %>%
  summary()
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)               -2.11824    0.02788 -75.971  < 2e-16 ***
#   median_trans               8.12840    0.08689  93.549  < 2e-16 ***
#   Relation_to_IslandN_Shelf  0.35235    0.02950  11.944  < 2e-16 ***
#   Relation_to_IslandN_Shore  0.33235    0.02007  16.561  < 2e-16 ***
#   Relation_to_IslandOpenSea  0.35733    0.01930  18.513  < 2e-16 ***
#   Relation_to_IslandS_Shelf  0.31476    0.03116  10.102  < 2e-16 ***
#   Relation_to_IslandS_Shore  0.35610    0.02126  16.746  < 2e-16 ***
#   GeneRelIntergenic         -0.16639    0.01493 -11.148  < 2e-16 ***
#   TssATRUE                   0.15682    0.01827   8.586  < 2e-16 ***
#   TssAFlnkTRUE               0.45341    0.01651  27.461  < 2e-16 ***
#   TxFlnkTRUE                 0.01230    0.02097   0.587  0.55753
# TxWkTRUE                   0.25240    0.01429  17.663  < 2e-16 ***
#   TxTRUE                    -0.01657    0.02075  -0.799  0.42453
# EnhGTRUE                   0.25193    0.02256  11.169  < 2e-16 ***
#   EnhTRUE                    0.32913    0.01399  23.528  < 2e-16 ***
#   ZNF.RptsTRUE               0.56330    0.03724  15.128  < 2e-16 ***
#   HetTRUE                   -0.07398    0.02513  -2.943  0.00325 **
#   TssBivTRUE                -0.10580    0.02520  -4.199 2.68e-05 ***
#   BivFlnkTRUE                0.32190    0.02201  14.624  < 2e-16 ***
#   EnhBivTRUE                -0.04924    0.02015  -2.444  0.01453 *
#   ReprPCTRUE                -0.02174    0.01627  -1.337  0.18132
# ReprPCWkTRUE               0.18674    0.01590  11.746  < 2e-16 ***
#   QuiesTRUE                  0.08442    0.01547   5.456 4.88e-08 ***
  

## Test variables explaining differences between CpG Types

## Mono_Inverse

runFullModel <- function(type) {
  as_tibble(methyAnnot) %>%
    mutate(CpG = Name) %>%
    right_join(CpGsSum, by = "CpG") %>%
    filter(Type != "Non-significant") %>%
    mutate(SelCpG = ifelse(Combined %in% type, "Sel", "Other"),
           SelCpG = factor(SelCpG, levels = c("Other", "Sel")),
           median_trans = -(median-0.5)^2) %>%
    glm(formula(paste("SelCpG ~ median_trans + Relation_to_Island + GeneRel +",
                      paste(chromStates, collapse = "+"))), 
        family = "binomial", .) %>%
    summary()
}
filt <- function(mod){
  mod$coefficient[mod$coefficient[, 4] < 0.05,]
}
runFullModel("Mono_Inverse")
filt(runFullModel("Mono_Inverse"))
runFullModel("Mono_Positive")
filt(runFullModel("Mono_Positive"))
runFullModel("Multi_Inverse")
filt(runFullModel("Multi_Inverse"))
runFullModel("Multi_Positive")
filt(runFullModel("Multi_Positive"))
runFullModel("Multi_Both")
filt(runFullModel("Multi_Both"))


runFullModel(c("Mono_Inverse", "Mono_Positive"))
runFullModel(c("Multi_Positive", "Mono_Positive"))

## Blood cell types ####
png("paper/CpGEnrichBloodTypes.png", width = 3000, height = 2000, res = 300)
FlowSorted.Blood.450k.compTable %>%
  mutate(CpG = rownames(.)) %>%
  select(CpG, Fstat) %>%
  right_join(CpGsSum, by = "CpG") %>% 
  as_tibble() %>%
  mutate(Combined = factor(Combined, 
                       levels = c("Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both", "Non-significant")
                       )
         ) %>%
  ggplot(aes(x = Combined, y = log10(Fstat), fill = Combined)) +
  geom_boxplot() +
  theme_bw() + 
  theme(legend.position = "none") +
  scale_x_discrete(name = "CpG Type") +
  scale_fill_manual(name = "CpG Type", values = c("#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73", "#FFFFFF")) 
dev.off()

FlowSorted.Blood.450k.compTable %>%
  mutate(CpG = rownames(.)) %>%
  select(CpG, Fstat) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Combined = factor(Combined, 
                           levels = c("Non-significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both")
                           )
  ) %>%
  lm(log10(Fstat) ~ Combined, .) %>%
  summary()

# Call:
#   lm(formula = log10(Fstat) ~ Combined, data = .)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -4.4132 -0.4053 -0.0312  0.3470  2.8148
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            0.628921   0.001050  598.77   <2e-16 ***
#   CombinedMono_Inverse   0.451689   0.005636   80.14   <2e-16 ***
#   CombinedMono_Positive  0.307423   0.006661   46.15   <2e-16 ***
#   CombinedMulti_Inverse  0.500182   0.008020   62.37   <2e-16 ***
#   CombinedMulti_Positive 0.336254   0.010538   31.91   <2e-16 ***
#   CombinedMulti_Both     0.471616   0.009946   47.42   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6212 on 384844 degrees of freedom
# (1568 observations deleted due to missingness)
# Multiple R-squared:  0.03674,   Adjusted R-squared:  0.03673
# F-statistic:  2936 on 5 and 384844 DF,  p-value: < 2.2e-16
# 


## Age variability ####
agedf <- read.csv2("data/DiffMethyAgeCpGs.csv", as.is = TRUE)
agedf <- agedf %>%
  as_tibble() %>%
  mutate(Dir = ifelse(as.numeric(beta8.avg) > as.numeric(beta0.avg), "Increasing", "Decreasing"), 
         CpG = ILMNID)

ageSum <- CpGsSum %>%
  dplyr::select(CpG, Direction) %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  left_join(dplyr::select(agedf, Dir, CpG), by = "CpG") %>%
  mutate(Dir = ifelse(is.na(Dir), "Constant", Dir)) %>%
  group_by(Direction, Dir) %>%
  summarize(n = n()) %>%
  spread(Dir, n) %>%
  mutate(Type0 = ifelse(Direction == "Non-significant", "Non-significant", "Significant"), 
         Changed = sum(Decreasing, Increasing), 
         tot = sum(Changed, Constant))

ageG <- c("Changed", "Decreasing", "Increasing")

allAge <- ageSum %>%
  group_by(Type0) %>%
  dplyr::select(-c(1:2)) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(ageG) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:3) %>%
  mutate(Type = "Significant")

typesAge <- lapply(c("Inverse", "Positive", "Both"), function(t){
  ageSum %>%
    filter(Direction %in% c(t, "Non-significant")) %>%
    group_by(Direction) %>%
    dplyr::select(-Type0) %>%
    summarize_all(list(sum)) %>%
   # arrange(Combined) %>%
    g2(ageG) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:length(ageG)) %>%
    mutate(Type = t)
})
combAge <- Reduce(rbind, typesAge)
combAge <- rbind(allAge, combAge)

png("paper/CpGEnrichAge.png", width = 3000, height = 2000, res = 300)
combAge %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Inverse", "Positive", "Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Change during chilhood", limits = ageG) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 
dev.off()

## Overlap with literature ####
#### EWAS Catalogue ####
ewasdf <- read.delim(gzfile("data/EWAS_Catalog_03-07-2019.txt.gz"))

## Restrict to results in blood and European
ewasdfF <- ewasdf %>%
  filter(grepl("Whole|Peripheral", Tissue)) %>%
  filter(!N_EUR == "-" & (N_EAS != "-" | N_SAS != "-" | N_AFR != "-" | N_AMR != "-" | N_OTH == "-"))

ewasCatSum <- ewasdfF %>%
  dplyr::select(CpG) %>%
  mutate(Present = "Present") %>%
  distinct() %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Present = ifelse(is.na(Present), "Absent", Present)) %>%
  as_tibble() 

ewasCatSum %>% 
  distinct() %>%
  summarize(all = sum(Present == "Present"),
            eQTMs = sum(Present == "Present" & Type != "Non-significant")) %>%
  mutate(eQTMs_cat = eQTMs/all)
# all eQTMs eQTMs_cat
# 143384 16083     0.112

t <- table(ewasCatSum$Present, ewasCatSum$Type != "Non-significant")
t[1]/t[2]/t[3]*t[4]


a <- ewasCatSum %>% 
  group_by(CpG) %>% 
  summarize(Cat = ifelse(n() == 2, "Exposure-Outcome", Category),
            Combined = unique(Combined)) %>%
  group_by(Cat, Combined) %>% 
  summarise(n = n()) %>%
  spread(key = Cat, value = n)
colSums(a[,-1])
# Exposure Exposure-Outcome             None          Outcome
# 142444              547           243034              393
colSums(a[,-1])/sum(colSums(a[, -1])) * 100
# Exposure Exposure-Outcome             None          Outcome
# 36.8626720        0.1415566       62.8940681        0.1017033


catalCpGs <- unique(ewasdfF$CpG)
CpGsSum %>%
  mutate(mQTL = CpG %in% catalCpGs) %>%
  group_by(Combined) %>%
  summarize(n = sum(mQTL),
            prop = mean(mQTL))

CpGsSum %>%
  mutate(mQTL = CpG %in% catalCpGs,
         Type = ifelse(Combined != "Non-significant", "Significant", "Non-significant")) %>%
  group_by(Type) %>%
  summarize(n = sum(mQTL),
            prop = mean(mQTL))

ewasCatSum %>% 
  group_by(CpG) %>% 
  summarize(Cat = ifelse(n() == 2, "Exposure-Outcome", Category),
            Combined = ifelse(unique(Combined) == "Non-significant", "Non-significant", "Significant")) %>%
  group_by(Cat, Combined) %>% 
  summarise(n = n()) %>%
  spread(key = Cat, value = n) %>%
  mutate(TotalCpGs = None + Exposure + `Exposure-Outcome` + Outcome, 
       catalogue = Exposure + `Exposure-Outcome` + Outcome,
       pCatalogue = catalogue/TotalCpGs,
       pExp = Exposure/TotalCpGs,
       pOutcome = Outcome/TotalCpGs,
       pComb = `Exposure-Outcome`/TotalCpGs)
  
ewasCatSum2 <-  ewasCatSum %>% 
  group_by(CpG) %>% 
  summarize(Cat = ifelse(n() == 2, "Exposure-Outcome", Category),
            Combined = unique(Combined)) %>%
  group_by(Combined, Cat) %>%
  summarize(n = n()) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant")) %>%
  spread(Cat, n) %>%
  mutate(tot = Exposure + `Exposure-Outcome` + Outcome + None) %>%
  select(-None)
  
catals <- c("Exposure", "Exposure-Outcome", "Outcome")


allCatal <- ewasCatSum2 %>%
  group_by(Type0) %>%
  select(-Combined) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(catals) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:3) %>%
  mutate(Type = "Significant")

typesCatal <- lapply(types, function(t){
  ewasCatSum2 %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    arrange(Combined) %>%
    g2(catals) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:3) %>%
    mutate(Type = t)
})
combCatal <- Reduce(rbind, typesCatal)
combCatal <- rbind(allCatal, combCatal)

png("paper/CpGEnrichCatalogue.png", width = 3000, height = 2000, res = 300)
combCatal %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "", limits = cats) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 
dev.off()

## Test specific exposures/outcomes ####
### Select categories with > 200 CpGs
expoCatal <- ewasdfF %>% 
  group_by(Exposure) %>%
  summarize(n = n()) %>%
  filter(n > 200) %>%
  pull(Exposure) %>%
  as.character()

## Exclude age and DNA methylation
expoCatal <- expoCatal[-grep("Age|me..ylation", expoCatal)]
## Exclude migration and disease
expoCatal <- expoCatal[!expoCatal %in% c("Migration in Italy", "Lupus nephritis vs lupus without nephritis")] 
## Group smoking in one category
expoCatal[expoCatal %in% c("Current versus never smoking", "Former versus never smoking")] <- "Smoking"
expoCatal[expoCatal == "Maternal smoking during pregnancy"] <- "Maternal smoking in pregnancy"
expoCatal <- unique(expoCatal)

dataCpGs <- CpGsSum$CpG
expoMat <- sapply(expoCatal, function(x){
  dataCpGs %in% ewasdfF[ewasdfF$Exposure == x, "CpG"]
})

## Test specific outcomes
### Select categories with > 200 CpGs
outCatal <- ewasdfF %>% 
  group_by(Outcome) %>%
  summarize(n = n()) %>%
  filter(n > 200) %>%
  pull(Outcome) %>%
  as.character()
## Exclude age and DNA methylation
outCatal <- outCatal[outCatal != "DNA methylation"]
outMat <- sapply(outCatal, function(x){
  dataCpGs %in% ewasdfF[ewasdfF$Outcome == x, "CpG"]
})
colnames(outMat)[1] <- "BMI_Outcome"

exps <- c(colnames(expoMat), colnames(outMat))

expoSum <- cbind(CpGsSum, expoMat, outMat) %>%
  as_tibble() %>%
  group_by(Combined) %>%
  select(-CpG, -Type, -Direction) %>%
  summarize_all(sum) %>%
  left_join(group_by(CpGsSum, Combined) %>% summarize(tot = n()), by = "Combined") %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"))

allExps <- expoSum %>%
  group_by(Type0) %>%
  select(-Combined) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(exps) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", seq_len(length(exps))) %>%
  mutate(Type = "Significant")

typesExps <- lapply(types, function(t){
  expoSum %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    arrange(Combined) %>%
    g2(exps) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value",seq_len(length(exps))) %>%
    mutate(Type = t)
})
combExps <- Reduce(rbind, typesExps)
combExps <- rbind(allExps, combExps)

png("paper/CpGEnrichCatalogueExposures.png", width = 3000, height = 2000, res = 300)
combExps %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both")),
         Group = ifelse(Region %in% c("Body mass index", "C-reactive protein", "Sex"), "Phenotype", 
                        ifelse(Region %in% c("Autoantibody production in systemic lupus erythematosus", "HIV infection", "Primary Sjogrens syndrome", "Rheumatoid arthritis", "Schizophrenia"), "Disease", 
                               ifelse(Region %in% c("BMI_Outcome", "Total serum IgE"), "Outcome", "Exposure")))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  facet_grid(~ Group, scales = "free", space = "free_x") +
  scale_x_discrete(name = "") +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() +
  theme(axis.text = element_text(angle = 45))
dev.off()

## EWAS Atlas #### 
## Download time: 27/11/2019
atAssoc <- read.delim("data/EWAS_Atlas_associations.tsv")
atStudy <- read.delim("data/EWAS_Atlas_studies.tsv")
atCohort <- read.delim("data/EWAS_Atlas_cohorts.tsv")




## Combined traits (web page does not allow to download all files)
# traitsL <- lapply(1:5, function(x) read.csv(paste0("data/tableExport", x, ".csv")))
# traits <- Reduce(rbind, traitsL)
traits <- read.csv("data/traitsTable.csv")

colnames(atStudy)[1] <- "Study.id"
colnames(atCohort)[c(2, 4, 11, 15)] <- c("Study.id", "Array", "Source", "Ancestry")
atlasDF <- atStudy %>%
  select(Study.id) %>%
  full_join(select(atCohort, Study.id, Array, Source, Ancestry), by = "Study.id") %>% 
  left_join(atAssoc, ., by ="Study.id") %>%
  as_tibble() %>%
  distinct()

## Select Europeans and blood and add trait type
atlasDFsel <- atlasDF %>%
  filter(Ancestry %in% c("Not reported", "European") &
           Source %in% c("whole blood", "peripheral blood")) %>%
  left_join(select(traits, Trait, Type), by = "Trait")


atlasSum <- CpGsSum %>%
  mutate(Changed = CpG %in% atlasDFsel$Probe.id) %>%
  distinct() 
atlasSum %>%
  summarize(all = sum(Changed),
            eQTMs = sum(Changed & Type != "Non-significant")) %>%
  mutate(eQTMs_cat = eQTMs/all)
# all eQTMs eQTMs_cat
# 54599  9547     0.175

t <- table(atlasSum$Changed, atlasSum$Type != "Non-significant")
t[1]/t[2]/t[3]*t[4]



## Get summary statistics
atlasSum %>% 
  group_by(Type0) %>%
  select(ends_with("in")) %>%
  summarize_all(sum)
  


allAtlas <- atlasSum %>%
  group_by(Type0) %>%
  select(ends_with("in"), ends_with("ou")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g(eff) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", seq_len(length(eff))) %>%
  mutate(Type = "Significant")

typesAtlas <- lapply(types, function(t){
  atlasSum %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(ends_with("in"), ends_with("ou")) %>%
    summarize_all(list(sum)) %>%
    g(eff) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", seq_len(length(eff))) %>%
    mutate(Type = t)
})
combAtlas <- Reduce(rbind, typesAtlas)
combAtlas <- rbind(allAtlas, combAtlas)


png("paper/CpGEnrichAtlas.png", width = 3000, height = 2000, res = 300)
combAtlas %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "EWAS atlas category", limits = eff) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 
dev.off()

## Test specific exposures/outcomes ####
### Select categories with > 200 CpGs
expoAtlas <- atlasDFsel %>% 
  group_by(Trait) %>%
  summarize(n = n()) %>%
  filter(n > 200) %>%
  pull(Trait) %>%
  as.character()

dataCpGs <- CpGsSum$CpG
expoAtlasMat <- sapply(expoAtlas, function(x){
  dataCpGs %in% atlasDFsel[atlasDFsel$Trait == x, "Probe.id", drop = TRUE]
})
expsA <- colnames(expoAtlasMat)
expoA <- cbind(CpGsSum, expoAtlasMat) %>%
  as_tibble() %>%
  group_by(Combined) %>%
  select(-CpG, -Type, -Direction) %>%
  summarize_all(sum) %>%
  left_join(group_by(CpGsSum, Combined) %>% summarize(tot = n()), by = "Combined") %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"))

allExpsAt <- expoA %>%
  group_by(Type0) %>%
  select(-Combined) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(expsA) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", seq_len(length(expsA))) %>%
  mutate(Type = "Significant")

typesExpsAt <- lapply(types, function(t){
  expoA %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    arrange(Combined) %>%
    g2(expsA) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value",seq_len(length(expsA))) %>%
    mutate(Type = t)
})
combExpsA <- Reduce(rbind, typesExpsAt)
combExpsA <- rbind(allExpsAt, combExpsA)
combExpsA <- combExpsA %>%
  spread(par, Value) 

expsSum <- combExpsA %>%
  mutate(Trait = Region) %>%
  group_by(Trait) %>%
  summarize(Sig = any(p.val < 1e-3)) %>%
  left_join(group_by(atlasDFsel, Trait) %>% summarize(n = n()), by = "Trait")

## Only plot exposures enriched in at least one CpG group
plotExps <- filter(expsSum, Sig) %>% pull(Trait)

a <- combExpsA %>%
  mutate(CpGType = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both")),
         Trait = Region) %>%
  select(-Type) %>%
  left_join(select(traits, Trait, Type), by = "Trait") %>%
  mutate(Group = ifelse(Type %in% c("cancer", "non-cancer disease"), "Disease", "Exposure")) %>%
  filter(Trait %in% plotExps) %>%
  ## Change OR of 0 to OR of 1 to avoid having meaningless bars in the plot
  mutate(OR = ifelse(OR == 0, 1, OR))
p1 <- a %>%
  filter(Group == "Disease") %>%
  ggplot(aes(x = Region, y = OR, fill = CpGType)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  facet_grid(~ Type, scales = "free", space = "free_x") +
  scale_x_discrete(name = "") +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() +
  theme(axis.text = element_text(angle = 45))
p2 <- a %>%
  filter(Group != "Disease") %>%
  ggplot(aes(x = Region, y = OR, fill = CpGType)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  facet_grid(~ Type, scales = "free", space = "free_x") +
  scale_x_discrete(name = "") +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() +
  theme(axis.text = element_text(angle = 45))

png("paper/CpGEnrichAtlasExposures.png", width = 3500, height = 4000, res = 300)
plot_grid(p1, p2, ncol = 1, labels = "")
dev.off()

