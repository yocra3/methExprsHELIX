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
library(dplyr)
library(cowplot)
library(tidyr)
library(hta20transcriptcluster.db)
library(FlowSorted.Blood.450k)
library(openxlsx)
library(GenomicRanges)

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
  mixed <- runTest(Data, algorithm = "weight01", statistic = "fisher")
  
  
  finTab <- GenTable(Data, 
                     w0 = mixed,
                     classic = class, 
                     orderBy = "w0",
                     ranksOf = "classic",
                     topNodes = length(score(class)))
  finTab <- subset(finTab, w0 < 0.01 & classic < 0.01)
  
  gos <- list(classic = class, mixed = mixed)
  list(gos = gos, table = finTab)
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


## Subtypes
subtypes <- c("Mono_Inverse", "Mono_Positive", "Multi_Both", "Multi_Inverse", 
              "Multi_Positive")
names(subtypes) <- subtypes
subtypes_go <- lapply(subtypes, runGO, annot = CpGsSum, df = df)

## Save as Rdata data and process in local
save(allGos, subtypes_go, file = "paper/GOobjects.Rdata")

## It does not work on server. Run locally.
library(GOfuncR)
library(dplyr)
server <- "//isg10174/data/WS_HELIX/HELIX_analyses/expr_met_SM/paper/"

load(paste0(server, "GOobjects.Rdata"))

addImmunityInfo <- function(tab){
  get_parent_nodes(tab$GO.ID) %>%
    group_by(child_go_id) %>%
    summarize(tag = ifelse(any(c("adaptive immune response", "lymphocyte activation", "antigen processing and presentation") %in% parent_name), "adaptive",
                           ifelse(any(c("innate immune response", "myeloid leukocyte activation", "myeloid leukocyte mediated immunity", "myeloid leukocyte differentiation", "myeloid leukocyte cytokine production") %in% parent_name), "innate",
                                  ifelse(any(c("response to cytokine", "immune system process") %in% parent_name), "general", "none")))) %>%
    mutate(immune = tag != "none", 
           GO.ID = child_go_id) %>%
    right_join(tab)  %>%
    select(GO.ID, Term, w0, classic, immune, tag)
}

allMod <- addImmunityInfo(allGos$tab)

subtypesTabs <- lapply(subtypes_go, function(x){
  addImmunityInfo(x$tab)
})
  
  
write.table(allMod[, c("GO.ID", "Term", "w0", "classic", "tag")], 
            file = paste0(server, "/GOsAllGenes.txt"), 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)

lapply(names(subtypesTabs), function(x){
  write.table(subtypesTabs[[x]][, c("GO.ID", "Term", "w0", "classic", "tag")], 
              file = paste0(server, "/GOs", x, ".txt"), 
              quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
  
})
table(allMod$tag)
a <- lapply(subtypesTabs, function(x){
  c(n = nrow(x), imm = sum(x$immune), immP = round(mean(x$immune)*100, 1), 
       adap = sum(x$tag == "adaptive"), innate = sum(x$tag == "innate"))
}) %>%
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

## Gene position ####
methGenePos <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, TSS200, TSS1500, UTR5, FirstExon, Body, UTR3, GeneRel) %>%
  right_join(CpGsSum, by = "CpG") %>%
  group_by(Combined) %>%
  summarize(Intergenicin = sum(GeneRel == "Intergenic"),
            Intergenicou = sum(GeneRel != "Intergenic"),
            TSS1500in = sum(TSS1500),
            TSS1500ou = sum(!TSS1500),
            TSS200in = sum(TSS200),
            TSS200ou = sum(!TSS200),
            UTR5in = sum(UTR5),
            UTR5ou = sum(!UTR5),
            FirstExonin = sum(FirstExon),
            FirstExonou = sum(!FirstExon),
            Bodyin = sum(Body),
            Bodyou = sum(!Body),
            UTR3in = sum(UTR3),
            UTR3ou = sum(!UTR3)) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"))

types <- c( "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both")
gpos <- c("Intergenic", "TSS1500", "TSS200", "UTR5", "FirstExon", "Body", "UTR3")

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

allPos <- methGenePos %>%
  group_by(Type0) %>%
  select(ends_with("in"), ends_with("ou")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g(gpos) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:7) %>%
  mutate(Type = "Significant")

typesPos <- lapply(types, function(t){
  methGenePos %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(ends_with("in"), ends_with("ou")) %>%
    summarize_all(list(sum)) %>%
    g(gpos) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:7) %>%
    mutate(Type = t)
})
combPos <- Reduce(rbind, typesPos)
combPos <- rbind(allPos, combPos)

png("paper/CpGEnrichGenePos.png", width = 3000, height = 2000, res = 300)
combPos %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Position with respect gene", 
                   limits = c("TSS1500", "TSS200", "UTR5", "FirstExon", "Body", "UTR3", "Intergenic")) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 
dev.off()


## Methylation levels ####
methMethLevels <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median_cat) %>%
  right_join(CpGsSum, by = "CpG") %>%
  group_by(Combined, median_cat) %>%
  summarize(n = n()) %>%
  spread(median_cat, n) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"), 
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

allMethLevs <- methMethLevels %>%
  group_by(Type0) %>%
  select(-Combined) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(cats) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:3) %>%
  mutate(Type = "Significant")

typesMethLevs <- lapply(types, function(t){
  methMethLevels %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    arrange(Combined) %>%
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
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both"))) %>%
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


## CpG Islands ####
methIsland <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, Relation_to_Island) %>%
  right_join(CpGsSum, by = "CpG") %>%
  group_by(Combined, Relation_to_Island) %>%
  summarize(n = n()) %>%
  spread(Relation_to_Island, n) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"), 
      tot = sum(Island, N_Shelf, N_Shore, OpenSea, S_Shelf, S_Shore))

islandStates <- c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "OpenSea")

allIsland <- methIsland %>%
  group_by(Type0) %>%
  select(-Combined) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(islandStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:length(islandStates)) %>%
  mutate(Type = "Significant")

typesIsland <- lapply(types, function(t){
  methIsland %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    arrange(Combined) %>%
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
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both"))) %>%
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

## Chromatin states ####
chromStates <- c("TssA", "TssAFlnk", "TxFlnk", "TxWk", "Tx", "EnhG", "Enh",
                 "ZNF.Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
                 "ReprPCWk", "Quies")
sum2 <- function(x) sum(!x)

methChromSt <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, eval(chromStates)) %>%
  right_join(CpGsSum, by = "CpG") %>%
  group_by(Combined) %>%
  summarize_at(chromStates, list(sum, sum2)) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"))


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
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
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
                       levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both")),
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
  facet_grid(~ Group, scales = "free", space = "free_x") +
  theme_bw() 
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

## Gemes Clusters ####
### PMID: 24656863
#### Usamos solo los del paper principal para tener más que los que están asociados a genotipos
cpgClus1 <- read.xlsx("data/CpGCluster_1.xlsx", startRow = 2)

cpgClus1GR <- makeGRangesFromDataFrame(cpgClus1, start.field = "Start.Position.(hg19)",
                                       end.field = "End.Position.(hg19)")

seqlevelsStyle(cpgClus1GR) <- "UCSC"
cpgClus <- methyAnnot %>%
  makeGRangesFromDataFrame(start.field = "pos", end.field = "pos") %>%
  findOverlaps(cpgClus1GR)
cpgClus <- rownames(methyAnnot)[from(cpgClus)]

genoClus <- read.xlsx("data/CpGClusterGeno.xlsx", startRow = 2)
genoClusvec <- unlist(strsplit(genoClus$CpG.Probe.Names, " "))

cegmeQTLs <- read.xlsx("data/meQTLs_CeGEMs.xlsx", startRow = 2)

CEGEMS <- CpGsSum %>%
  mutate(clus = CpG %in% cpgClus,
         meql = CpG %in% unique(cegmeQTLs$CpG),
    geno = CpG %in% genoClusvec) %>%
  group_by(Combined) %>%
  summarize(clusin = sum(clus),
            clusou = sum2(clus),
            meQTLin = sum(meql),
            meQTLou = sum2(meql),
            genoin = sum(geno),
            genoou = sum2(geno)) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-Significant", "Significant"))
  
cegs <- c("clus", "meQTL", "geno")

allCeg <- CEGEMS %>%
  group_by(Type0) %>%
  select(ends_with("in"), ends_with("ou")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g(cegs) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:3) %>%
  mutate(Type = "Significant")

typesCeg <- lapply(types, function(t){
  CEGEMS %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(ends_with("in"), ends_with("ou")) %>%
    summarize_all(list(sum)) %>%
    g(cegs) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:3) %>%
    mutate(Type = t)
})
combCeg <- Reduce(rbind, typesCeg)
combCeg <- rbind(allCeg, combCeg)

png("paper/CpGEnrichCEGEMS.png", width = 3000, height = 2000, res = 300)
combCeg %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "GeMes classification", 
                   limits = cegs,
                   labels = c("Contiguous Cluster", "meQTL", "Genotype-Controlled Cluster")) +
  scale_fill_manual(name = "CpG Type", values = c("#999999", "#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73")) +
  theme_bw() 
dev.off()

## Age variability ####
agedf <- read.csv2("data/DiffMethyAgeCpGs.csv", as.is = TRUE)
agedf <- agedf %>%
  as_tibble() %>%
  mutate(Dir = ifelse(as.numeric(beta8.avg) > as.numeric(beta0.avg), "Increasing", "Decreasing"), 
         CpG = ILMNID)

ageSum <- CpGsSum %>%
  select(CpG, Combined) %>%
  left_join(select(agedf, Dir, CpG), by = "CpG") %>%
  mutate(Dir = ifelse(is.na(Dir), "Constant", Dir)) %>%
  group_by(Combined, Dir) %>%
  summarize(n = n()) %>%
  spread(Dir, n) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"), 
         Changed = sum(Decreasing, Increasing), 
         tot = sum(Changed, Constant))

ageG <- c("Changed", "Decreasing", "Increasing")

allAge <- ageSum %>%
  group_by(Type0) %>%
  select(-c(1:2)) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g2(ageG) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:3) %>%
  mutate(Type = "Significant")

typesAge <- lapply(types, function(t){
  ageSum %>%
    filter(Combined %in% c(t, "Non-significant")) %>%
    group_by(Combined) %>%
    select(-Type0) %>%
    summarize_all(list(sum)) %>%
    arrange(Combined) %>%
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
  mutate(Type = factor(Type, levels = c("Significant", "Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both"))) %>%
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
ewasdfF <- mutate(ewasdf, Category = ifelse(Outcome %in% c("DNA methylation", "Changes in DNA methylation"), "Exposure", "Outcome")) %>%
  filter(grepl("Whole|Peripheral", Tissue)) %>%
  filter(!N_EUR == "-" & (N_EAS != "-" | N_SAS != "-" | N_AFR != "-" | N_AMR != "-" | N_OTH == "-"))

ewasCatSum <- ewasdfF %>%
  select(CpG, Category) %>%
  mutate(Present = "Present") %>%
  distinct() %>%
  right_join(CpGsSum, by = "CpG") %>%
  mutate(Category = ifelse(is.na(Category), "None", Category), 
         Present = ifelse(is.na(Present), "Absent", Present)) %>%
  as_tibble()

ewasCatSum %>% 
  select(-Category) %>%
  distinct() %>%
  summarize(n = sum(Present == "Present"),
            mean = mean(Present == "Present"))
# n  
# <int> <dbl>
# 143384 0.371
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
## Download time: 06/09/2019
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
  mutate(Changed = CpG %in% atlasDFsel$Probe.id,
         Behavior = CpG %in% filter(atlasDFsel, Type == "behavior")$Probe.id, 
         Cancer = CpG %in% filter(atlasDFsel, Type == "cancer")$Probe.id,
         Environment = CpG %in% filter(atlasDFsel, Type == "environmental factor")$Probe.id, 
         Disease = CpG %in% filter(atlasDFsel, Type == "non-cancer disease")$Probe.id, 
         Phenotype = CpG %in% filter(atlasDFsel, Type == "phenotype")$Probe.id) %>%
  group_by(Combined) %>%
  summarize(Changedin = sum(Changed),
            Changedou = sum(!Changed),
            Cancerin = sum(Cancer),
            Cancerou = sum(!Cancer),
            Environmentin = sum(Environment),
            Environmentou = sum(!Environment),
            Behaviorin = sum(Behavior),
            Behaviorou = sum(!Behavior),
            Diseasein = sum(Disease),
            Diseaseou = sum(!Disease),
            Phenotypein = sum(Phenotype),
            Phenotypeou = sum(!Phenotype)) %>%
  mutate(Type0 = ifelse(Combined == "Non-significant", "Non-significant", "Significant"))
eff <- c("Changed", "Behavior", "Disease", "Environment", "Cancer", "Phenotype")

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

