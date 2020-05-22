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
server <- "/run/user/1000/gvfs/sftp:host=isgws06.isglobal.org,port=2222,user=cruiz/home/isglobal.lan/cruiz/data/WS_HELIX/HELIX_analyses/expr_met_SM/paper/"

load(paste0(server, "GOobjects.Rdata"))
load(paste0(server, "snpGOs.Rdata"))

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
    summarize(parent = sumF(parent_name, topGOs = top_GOs), 
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
snpMod <- addImmunityInfo(subset(snpGos$table, as.numeric(w0) < 0.001))

subtypesTabs <- lapply(subtypes_go, function(x){
  subset(addImmunityInfo(x$tab), as.numeric(w0) < 0.001)
})
  
  
write.table(allMod[, c("GO.ID", "GO_term", "w0", "parent", "immune")], 
            file = paste0(server, "/GOsAllGenes.txt"), 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
write.table(snpMod[, c("GO.ID", "GO_term", "w0", "parent", "immune")], 
            file = paste0(server, "/GOsGenesSNPs.txt"), 
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
types <- c("eQTM", "Positive", "Inverse")

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
fsum <- function(x){
  if ( is.factor(x))
    "eQTM"
  else {
    sum(x)
}
    }

methMethLevels <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median_cat) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  group_by(Direction, median_cat) %>%
  summarize(n = n()) %>%
  spread(median_cat, n) %>%
  ungroup() %>%
  rbind(filter(., Direction != "Non-significant") %>% 
            summarize_all(fsum)) %>%
  mutate(Direction = as.character(Direction),
         Direction = ifelse(is.na(Direction), "eQTM", Direction),
         tot = low + medium + high)

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

methLevs_prop_plot <- methMethLevels %>%
  filter(Direction != "Both") %>%
  mutate(Low = low/tot,
         Medium = medium/tot,
         High = high/tot) %>%
  select(Direction, Low, Medium, High) %>%
  gather(categories, proportion, 2:4) %>%
  ungroup() %>%
  mutate(Direction = as.character(Direction),
         Direction = ifelse(Direction == "Non-significant", "non-eQTM", Direction),
        Direction = factor(Direction,  levels = c("non-eQTM", "eQTM", "Inverse", "Positive")),
        categories = factor(categories, levels = c("Low", "Medium", "High"))) %>%
  ggplot(aes(x = categories, y = proportion*100, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(name = "Methylation Levels") +
  scale_y_continuous(name = "Proportion of CpGs (%)") +
  scale_fill_manual(name = "CpG Type", values = c("#000000", "#999999", "#E69F00", "#009E73")) +
  theme_bw()

png("paper/CpGFreqsMethLevels.png", width = 3000, height = 2000, res = 300)
methLevs_prop_plot
dev.off()

typesMethLevs <- lapply(types, function(t){
  rbind(filter(methMethLevels, Direction == t),
        filter(methMethLevels, Direction == "Non-significant")) %>%
    g2(cats) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:3) %>%
    mutate(Type = t)
})
combMethLevs <- Reduce(rbind, typesMethLevs)

methLevsPlot <- combMethLevs %>%
  spread(par, Value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Methylation Levels", limits = cats, labels = tools::toTitleCase(cats)) +
  scale_fill_manual(name = "eQTM Type", values = c("#999999", "#E69F00", "#009E73")) +
  theme_bw() 

png("paper/CpGEnrichMethLevels.png", width = 3000, height = 2000, res = 300)
methLevsPlot
dev.off()


## CpG Islands ####
methIsland <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, Relation_to_Island) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  group_by(Direction, Relation_to_Island) %>%
  summarize(n = n()) %>%
  spread(Relation_to_Island, n) %>%
  ungroup() %>%
  rbind(filter(., Direction != "Non-significant") %>% 
          summarize_all(fsum)) %>%
  mutate(Direction = as.character(Direction),
         Direction = ifelse(is.na(Direction), "eQTM", Direction),
         tot = Island + N_Shelf + N_Shore + OpenSea + S_Shelf + S_Shore)

  
islandStates <- c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "OpenSea")


Island_prop_plot <- methIsland %>%
  filter(Direction != "Both") %>%
  mutate(N_Shelf = N_Shelf/tot,
         N_Shore = N_Shore/tot,
         Island = Island/tot,
         S_Shore = S_Shore/tot,
         S_Shelf = S_Shelf/tot,
         OpenSea = OpenSea/tot) %>%
  select(-tot) %>%
  gather(categories, proportion, 2:(2+length(islandStates) - 1)) %>%
  ungroup() %>%
  mutate(Direction = as.character(Direction),
         Direction = ifelse(Direction == "Non-significant", "non-eQTM", Direction),
         Direction = factor(Direction,  levels = c("non-eQTM", "eQTM", "Inverse", "Positive")),
         categories = factor(categories, levels = islandStates)) %>%
  ggplot(aes(x = categories, y = proportion*100, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(name = "Methylation Levels") +
  scale_y_continuous(name = "Proportion of CpGs (%)") +
  scale_fill_manual(name = "CpG Type", values = c("#000000", "#999999", "#E69F00", "#009E73")) +
  theme_bw()

png("paper/CpGFreqsIsland.png", width = 3000, height = 2000, res = 300)
Island_prop_plot
dev.off()

typesIsland <- lapply(types, function(t){
    rbind(filter(methIsland, Direction == t),
          filter(methIsland, Direction == "Non-significant")) %>%
    g2(islandStates) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:length(islandStates)) %>%
    mutate(Type = t)
})
combIsland <- Reduce(rbind, typesIsland)

cpgPosPlot <- combIsland %>%
  spread(par, Value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "CpG Island position", limits = islandStates) +
  scale_fill_manual(name = "eQTM Type", values = c("#999999", "#E69F00", "#009E73")) +
  theme_bw() 

png("paper/CpGEnrichIsland.png", width = 3000, height = 2000, res = 300)
cpgPosPlot
dev.off()

# CpG Islands vs Methylation levels
isl_all_meth <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, Relation_to_Island, median) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Significant = ifelse(Type == "Non-significant", "Non-eQTM", "eQTM")) %>%
  ggplot(aes(x = Relation_to_Island, y = median, fill = Significant)) +
  geom_violin() +
  scale_x_discrete(name = "CpG island", limits = islandStates) +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "", values = c("#999999", "#FFFFFF")) +
  theme_bw()

png("paper/Islands_methLevs_violin.png", width = 3000, height = 1000, res = 300)
isl_all_meth
dev.off()


## Chromatin states ####
chromStates <- c("TssA", "TssAFlnk", "TxFlnk", "TxWk", "Tx", "EnhG", "Enh",
                 "ZNF.Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
                 "ReprPCWk", "Quies")
sum2 <- function(x) sum(!x)


methChromStIni <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, eval(chromStates)) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = ifelse(Type == "Non-significant", "non-eQTM", "eQTM")) %>%
  group_by(Direction) %>%
  summarize_at(chromStates, list(sum = sum, sum2 = sum2, prop = mean))

methChromSt <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, eval(chromStates)) %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  group_by(Direction) %>%
  summarize_at(chromStates, list(sum = sum, sum2 = sum2, prop = mean))  %>%
  filter(Direction != "Non-significant") %>%
  rbind(methChromStIni)


chromSt_prop_plot <- methChromSt %>%
  filter(Direction != "Both") %>%
  select(Direction, ends_with("prop")) %>%
  gather(categories, proportion, 2:(2+length(chromStates) - 1)) %>%
  ungroup() %>%
  mutate(Direction = factor(Direction,  levels = c("non-eQTM", "eQTM", "Inverse", "Positive")),
         categories = gsub("_prop", "", categories),
         Group = factor(ifelse(categories %in% c("TssA", "TssAFlnk"), "TssProxProm",
                               ifelse(categories %in% c("Tx", "TxWk"), "ActTrans", 
                                      ifelse(categories %in% c("Enh", "EnhG"), "Enhancer", 
                                             ifelse(categories %in% c("TssBiv", "BivFlnk", "EnhBiv"), "BivReg", 
                                                    ifelse(categories %in% c("ReprPC", "ReprPCWk"), "ReprPoly", categories)
                                             )
                                      )
                               )
         ), 
         levels = c("TssProxProm", "TxFlnk", "ActTrans", "Enhancer", "ZNF.Rpts", "BivFlnk", "BivReg", "Het", "ReprPoly", "Quies")
         ),
         categories = factor(categories, levels = chromStates)) %>%
  ggplot(aes(x = categories, y = proportion*100, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  facet_wrap(~ Group, scales = "free_x") +  
  scale_y_continuous(name = "Proportion of CpGs (%)") +
  scale_fill_manual(name = "CpG Type", values = c("#000000", "#999999", "#E69F00", "#009E73")) +
  theme_bw()



typesChromSt <- lapply(types, function(t){
  rbind(filter(methChromSt, Direction == t),
        filter(methChromSt, Direction == "non-eQTM")) %>%
    select(ends_with("sum"), ends_with("sum2")) %>%
    g(chromStates, cols = c("_sum", "_sum2")) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:15) %>%
    mutate(Type = t)
})
combChromSt <- Reduce(rbind, typesChromSt)

png("paper/CpGEnrichChromStates.png", width = 3000, height = 2000, res = 300)
chromStatesPlot <- combChromSt %>%
  spread(par, Value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive")),
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
  scale_fill_manual(name = "eQTM Type", values = c("#999999", "#E69F00", "#009E73")) +
  facet_wrap(~ Group, scales = "free_x") +
  theme_bw() 
chromStatesPlot
dev.off()

## Chromatin states vs methylation levels
chrom_all_meth <-   as_tibble(methyAnnot) %>%
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
  geom_violin() +
  facet_wrap(~ Group, scales = "free_x") +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "", values = c("#999999", "#FFFFFF")) +
  theme_bw()
png("paper/medianMethvsChromState.png", width = 3000, height = 3000, res = 300)
chrom_all_meth
dev.off()

### Combined plot
png("paper/enrich_combined.png", width = 5000, height = 3000, res = 300)
plot_grid(
  plot_grid(cpgPosPlot, methLevsPlot, labels = c("A", "C"), nrow = 1),
  chromStatesPlot, 
  labels = c("", "B"), ncol = 1, rel_heights = c(1, 2)
)
dev.off()

png("paper/featuresProp_combined.png", width = 5000, height = 3000, res = 300)
plot_grid(
  plot_grid(Island_prop_plot, methLevs_prop_plot, labels = c("A", "C"), nrow = 1),
  chromSt_prop_plot, 
  labels = c("", "B"), ncol = 1, rel_heights = c(1, 2)
)
dev.off()

png("paper/enrich_methLevs_violin.png", width = 5000, height = 3000, res = 300)
plot_grid(isl_all_meth, chrom_all_meth, labels = LETTERS[1:2], ncol = 1, rel_heights = c(1, 2))
dev.off()



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





## Overlap with literature ####
## Age variability ####
agedf <- read.csv2("data/DiffMethyAgeCpGs.csv", as.is = TRUE)
agedf <- agedf %>%
  as_tibble() %>%
  mutate(Dir = ifelse(as.numeric(beta8.avg) > as.numeric(beta0.avg), "Increased", "Decreased"), 
         CpG = ILMNID)

ageSum <- CpGsSum %>%
  dplyr::select(CpG, Direction) %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  left_join(dplyr::select(agedf, Dir, CpG), by = "CpG") %>%
  mutate(Dir = ifelse(is.na(Dir), "Constant", Dir)) %>%
  group_by(Direction, Dir) %>%
  summarize(n = n()) %>%
  spread(Dir, n) %>%
  ungroup() %>%
  rbind(filter(., Direction != "Non-significant") %>% 
          summarize_all(fsum)) %>%
  mutate(Direction = as.character(Direction),
         Direction = ifelse(is.na(Direction), "eQTM", Direction),
         Variable = Decreased + Increased,
         tot = Variable + Constant)

ageG <- c("Variable", "Decreased", "Increased")


typesAge <- lapply(types, function(t){
  rbind(filter(ageSum, Direction == t),
        filter(ageSum, Direction == "Non-significant")) %>%
    g2(ageG) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:3) %>%
    mutate(Type = t)
})
combAge <- Reduce(rbind, typesAge)

png("paper/CpGEnrichAge.png", width = 3000, height = 2000, res = 300)
age_var <- combAge %>%
  spread(par, Value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Methylation during childhood", limits = ageG) +
  scale_fill_manual(name = "eQTM Type", values = c("#999999", "#E69F00", "#009E73")) +
  theme_bw() 
age_var
dev.off()

tmp <- CpGsSum %>%
  dplyr::select(CpG, Direction) %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  left_join(dplyr::select(agedf, Dir, CpG), by = "CpG") %>%
  mutate(Dir = ifelse(is.na(Dir), "Constant", Dir)) 

png("paper/CpGEnrichAge_methLevels.png", width = 3000, height = 2000, res = 300)
age_meth <-  as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median) %>%
  right_join(tmp, by = "CpG") %>%
  mutate(Significant = ifelse(Direction == "Non-significant", "Non-eQTM", "eQTM")) %>%
  ggplot(aes(x = Dir, y = median, fill = Significant)) +
  geom_violin() +
  scale_x_discrete(name = "Methylation during childhood") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "", values = c("#999999", "#FFFFFF")) +
  theme_bw()
age_meth
dev.off()

#### EWAS Catalog ####
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

ewasCatSum2 <-  ewasCatSum %>% 
  group_by(CpG) %>% 
  group_by(Present, Direction) %>%
  summarize(n = n()) %>%
  spread(Present, n) %>%
  mutate(Direction = factor(Direction)) %>%
  rbind(filter(., Direction != "Non-significant") %>% 
          summarize_all(fsum))  %>%
  mutate(Direction = as.character(Direction)) %>%
  select(Direction, Present, Absent)

typesCatal <- lapply(types, function(t){
  rbind(filter(ewasCatSum2, Direction == t),
        filter(ewasCatSum2, Direction == "Non-significant")) %>%
    getOR(cols = 2:3) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    mutate(Type = t, 
           dataset = "EWAS Catalog" )
})
combCatal <- Reduce(rbind, typesCatal)


png("paper/CpGEnrichCatalogue.png", width = 3000, height = 2000, res = 300)
combCatal %>%
  spread(par, value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Type, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "") +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#E69F00", "#009E73")) +
  theme_bw() 
dev.off()

catal_meth <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median) %>%
  right_join(ewasCatSum, by = "CpG") %>%
  mutate(Significant = ifelse(Type == "Non-significant", "Non-eQTM", "eQTM")) %>%
  ggplot(aes(x = Present, y = median, fill = Significant)) +
  geom_violin() +
  scale_x_discrete(name = "Presence in EWAS catalogue") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "", values = c("#999999", "#FFFFFF")) +
  theme_bw()



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
  group_by(Type) %>%
  select(ends_with("in")) %>%
  summarize_all(sum)
  

atlasSum2 <-  atlasSum %>% 
  mutate(Present = ifelse(Changed, "Present", "Absent")) %>%
  group_by(Present, Direction) %>%
  summarize(n = n()) %>%
  spread(Present, n) %>%
  mutate(Direction = factor(Direction)) %>%
  rbind(filter(., Direction != "Non-significant") %>% 
          summarize_all(fsum))  %>%
  mutate(Direction = as.character(Direction)) %>%
  select(Direction, Present, Absent)

  
typesAtlas <- lapply(types, function(t){
  rbind(filter(atlasSum2, Direction == t),
        filter(atlasSum2, Direction == "Non-significant")) %>%
    getOR(cols = 2:3) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    mutate(Type = t, 
           dataset = "EWAS Atlas")
})
combAtlas <- Reduce(rbind, typesAtlas)


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

png("paper/CpGEnrichEWASdbs.png", width = 3000, height = 2000, res = 300)
ewas_db <- combCatal %>%
  rbind(combAtlas) %>%
  spread(par, value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Type, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "") +
  scale_fill_manual(name = "eQTM Type", values = c("#999999", "#E69F00", "#009E73")) +
  theme_bw() +
  facet_grid(~  dataset)
ewas_db
dev.off()

png("paper/CpGEnrichEWASdbs_methLevels.png", width = 3000, height = 2000, res = 300)
as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median) %>%
  right_join(ewasCatSum, by = "CpG") %>%
  mutate(Significant = ifelse(Type == "Non-significant", "Non-eQTM", "eQTM"),
         Catalogue = Present) %>%
  select(Catalogue, Significant, CpG, median) %>%
  right_join(atlasSum, by = "CpG") %>%
  mutate(Atlas = ifelse(Changed, "Present", "Absent")) %>%
  select(Catalogue, Atlas, Significant, median) %>%
  gather(Database, Status, 1:2) %>%
  mutate(Database = ifelse(Database == "Catalogue", "EWAS Catalog", "EWAS Atlas")) %>%
  ggplot(aes(x = Status, y = median, fill = Significant)) +
  geom_violin() +
  scale_x_discrete(name = "Presence in catalogue") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "", values = c("#999999", "#FFFFFF")) +
  facet_grid(~ Database) +
  theme_bw()
dev.off()

### Make combined plot with age variant CpGs
png("paper/CpGEnrich_EWASdbs_age.png", width = 3000, height = 2500, res = 300)
plot_grid(ewas_db, age_var, ncol = 1, labels = "AUTO")
dev.off()