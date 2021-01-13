###############################################################################
# Code for used plots
###############################################################################

# Models comparison
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(UpSetR)

load("results/preprocessFiles/gexpAnnotation.Rdata")

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df
featsU <- featStatsDF

load("results/MethComBatExpResidualsCellAdj/allres_simP_cpgs.Rdata")
modC <- df
featsC <- featStatsDF

## CpGs identified
sigFU <- subset(featsU, p.val.adj < 0.05)$feat
sigFC <- subset(featsC, p.val.adj < 0.05)$feat
sigCom <- intersect(sigFU, sigFC)

png("compareModelsCpG.png")
upset(fromList(list("Model A" = sigFU, "Model B" = sigFC)), order.by = "freq", 
      set_size.angles = 45, text.scale = 3.5)
dev.off()

# Multiple testing comparison
library(dplyr)
library(UpSetR)

load("results/preprocessFiles/gexpAnnotation.Rdata")
expAnnot$TC <- expAnnot$transcript_cluster_id

## No cell
load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df %>%
  as_tibble() %>%
  mutate(adj.p.value = p.adjust(p.value, "BH"),
         sigPair_BH = adj.p.value < 0.05,
         adj.p.value.BF = p.adjust(p.value, "bonferroni"),
         sigPair_BF = adj.p.value.BF < 0.05,
         pair = paste(CpG, TC))
featsU <- featStatsDF

## CpGs identified
permsCp <- unique(modU[modU$sigPair, ]$CpG)
BHCp <- unique(modU[modU$sigPair_BH, ]$CpG)
BFCp <- unique(modU[modU$sigPair_BF, ]$CpG)

png("CpGsComparisonPval.png", height = 600)
upset(fromList(list(Permutation = permsCp, BH = BHCp, Bonferroni = BFCp)), 
      order.by = "freq", set_size.angles = 45, text.scale = 3.5)
dev.off()


## Not cell adjusted model ####
library(ggplot2)
library(dplyr)
library(tidyr)

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/gexpAnnotation.Rdata")


df_comp <- inner_join(df, overDF, by = c("CpG", "TC"))
allCpGs <- unique(df$CpG)
sigdf <- df[df$sigPair, ]
sigCpGs <- unique(sigdf$CpG)
CpGsSum <- df %>%
  group_by(CpG) %>%
  summarise(Type = ifelse(sum(sigPair) == 0, "Non-significant",
                          ifelse(sum(sigPair) == 1, "Mono", "Multi")),
            Direction = ifelse(sum(sigPair) == 0, "Non-significant",
                               ifelse(all(FC[sigPair] > 0), "Positive", 
                                      ifelse(all(FC[sigPair] < 0), "Inverse", "Both"))))


## Distance distribution CpG-TC pairs ####
png("PairDistr.png", width = 22.68, height = 13.45, units = "cm", res = 300)
ggplot(df_comp, aes(x = Distance, color = sigPair)) + geom_density() + 
  theme_bw() + 
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "Density", breaks = NULL) +
  scale_color_discrete(name = "", labels = c("Non-significant", "Significant")) +
  theme(plot.title = element_text(family = "Calibri", face = "bold", size = 35),
        axis.title = element_text(family = "Calibri", size = 30), 
        axis.text = element_text(family = "Calibri", size = 25),
        legend.title = element_text(family = "Calibri", size = 30), 
        legend.text = element_text(family = "Calibri", size = 20), 
        strip.text.x = element_text(family = "Calibri", size = 30))
dev.off()

## Enrichments ####
### Gene position ####
rownames(methyAnnot) <- methyAnnot$Name
methyAnnot <- methyAnnot[allCpGs, ]
methyAnnot$assoc <- methyAnnot$Name %in% sigCpGs
methyAnnot$GeneRel <- ifelse(methyAnnot$UCSC_RefGene_Name == "", "Intergenic", "Genic")
methyAnnot <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  left_join(CpGsSum, by = "CpG")

## Create summary tibble
posSum <- methyAnnot %>%
  group_by(Type, Direction) %>%
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
  mutate(Type0 = ifelse(Type == "Non-significant", "Non-significant", "Significant"))

## Define vars and function
gpos <- c("Intergenic", "TSS1500", "TSS200", "UTR5", "FirstExon", "Body", "UTR3")
getOR <- function(cols, df){
  t <- data.matrix(df[, cols])
  or <- t[1, 1]/(t[1, 2])/(t[2, 1])*t[2, 2]
  p.val <- chisq.test(t)$p.value
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
}

g <- function(x, gpos){
  sapply(gpos, function(y) getOR(paste0(y, c("in", "ou")), df = x))
}

sig <- posSum %>%
  group_by(Type0) %>%
  select(ends_with("in"), ends_with("ou")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g(gpos) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:7) %>%
  mutate(Type = "Significant")

neg <- posSum %>%
  filter(Direction %in% c("Inverse", "Non-significant")) %>%
  group_by(Direction) %>%
  select(ends_with("in"), ends_with("ou")) %>%
  summarize_all(list(sum)) %>%
  g(gpos) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:7) %>%
  mutate(Type = "Inverse")

pos <- posSum %>%
  filter(Direction %in% c("Positive", "Non-significant")) %>%
  group_by(Direction) %>%
  select(ends_with("in"), ends_with("ou")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Direction)) %>%
  g(gpos) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:7) %>%
  mutate(Type = "Positive")

both <- posSum %>%
  filter(Direction %in% c("Both", "Non-significant")) %>%
  group_by(Direction) %>%
  select(ends_with("in"), ends_with("ou")) %>%
  summarize_all(list(sum)) %>%
  g(gpos) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:7) %>%
  mutate(Type = "Both")

png("EnrichGenePos.png", width = 24.71, height = 12.25, units = "cm", res = 300)
rbind(sig, neg, pos, both) %>%
  spread(par, Value) %>%
  mutate(Type = factor(Type, levels = c("Significant", "Inverse", "Positive", "Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "", 
                   limits = c("TSS1500", "TSS200", "UTR5", "FirstExon", "Body", "UTR3", "Intergenic")) +
  theme_bw() +
  theme(plot.title = element_text(family = "Calibri", face = "bold", size = 35),
        axis.title = element_text(family = "Calibri", size = 25), 
        axis.text = element_text(family = "Calibri", size = 15),
        legend.title = element_text(family = "Calibri", size = 25), 
        legend.text = element_text(family = "Calibri", size = 15), 
        strip.text.x = element_text(family = "Calibri", size = 25))
dev.off()


## CpG Island ####
islandStates <- c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "OpenSea")

islSum <- methyAnnot %>%
  group_by(Type, Direction, Relation_to_Island) %>%
  summarize(n = n()) %>%
  spread(Relation_to_Island, n) %>%
  mutate(Type0 = ifelse(Type == "Non-significant", "Non-significant", "Significant"), 
         tot = sum(Island, N_Shelf, N_Shore, OpenSea, S_Shelf, S_Shore))

## Define vars and function
getOR <- function(col, df){
  t <- data.matrix(df[, c(col, colnames(df)[ncol(df)])])
  t[, 2] <- t[, 2] - t[, 1]
  or <- t[1, 1]/(t[1, 2])/(t[2, 1])*t[2, 2]
  p.val <- chisq.test(t)$p.value
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
}

g <- function(x, gpos){
  sapply(gpos, function(y) getOR(y, df = x))
}


sigI <- islSum %>%
  group_by(Type0) %>%
  select(-c(1:2)) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g(islandStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:6) %>%
  mutate(Type = "Significant")

NegI <- islSum %>%
  filter(Direction %in% c("Inverse", "Non-significant")) %>%
  group_by(Direction) %>%
  select(-starts_with("Type")) %>%
  summarize_all(list(sum))  %>%
  g(islandStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:6) %>%
  mutate(Type = "Inverse")

PosI <- islSum %>%
  filter(Direction %in% c("Positive", "Non-significant")) %>%
  group_by(Direction) %>%
  select(-starts_with("Type")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Direction)) %>%
  g(islandStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:6) %>%
  mutate(Type = "Positive")

BothI <- islSum %>%
  filter(Direction %in% c("Both", "Non-significant")) %>%
  group_by(Direction) %>%
  select(-starts_with("Type")) %>%
  summarize_all(list(sum)) %>%
  g(islandStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:6) %>%
  mutate(Type = "Both")


png("EnrichIsland.png", width = 24.24, height = 12.02, units = "cm", res = 300)
rbind(sigI, NegI, PosI, BothI) %>%
  spread(par, Value) %>%
  mutate(p.val.thres = ifelse(p.val > 0.05, "P > 0.05", 
                              ifelse(p.val < 1e-3, "P < 0.001", "P < 0.05")),
         Type = factor(Type, levels = c("Significant", "Inverse", "Positive", "Both"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  scale_fill_discrete(name = "") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "CpG Island", limits  = islandStates) +
  theme_bw() +
    theme(plot.title = element_text(family = "Calibri", face = "bold", size = 35),
          axis.title = element_text(family = "Calibri", size = 25), 
          axis.text = element_text(family = "Calibri", size = 15),
          legend.title = element_text(family = "Calibri", size = 25), 
          legend.text = element_text(family = "Calibri", size = 15), 
          strip.text.x = element_text(family = "Calibri", size = 25))
dev.off()

## Chromatin states ####
sum2 <- function(x) sum(!x)

chromStates <- c("TssA", "TssAFlnk", "TxFlnk", "TxWk", "Tx", "EnhG", "Enh",
                 "ZNF.Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
                 "ReprPCWk", "Quies")

chromSum <- methyAnnot %>%
  group_by(Type, Direction) %>%
  summarize_at(chromStates, funs(sum, sum2)) %>%
  mutate(Type0 = ifelse(Type == "Non-significant", "Non-significant", "Significant"))

## Define vars and function
getOR <- function(cols, df){
  t <- data.matrix(df[, cols])
  or <- t[1, 1]/(t[1, 2])/(t[2, 1])*t[2, 2]
  p.val <- chisq.test(t)$p.value
  ORl <- log(or) 
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])  
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl),      ORM = exp(ORl + 1.96*SEl))
}

g <- function(x, gpos){
  sapply(gpos, function(y) getOR(paste0(y, c("_sum", "_sum2")), x))
}

sigC <- chromSum %>%
  group_by(Type0) %>%
  select(ends_with("sum"), ends_with("sum2")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Type0)) %>%
  g(chromStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:15) %>%
  mutate(Type = "Significant")

negC <- chromSum %>%
  filter(Direction %in% c("Inverse", "Non-significant")) %>%
  group_by(Direction) %>%
  select(ends_with("sum"), ends_with("sum2")) %>%
  summarize_all(list(sum)) %>%
  g(chromStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:15) %>%
  mutate(Type = "Inverse")

posC <- chromSum %>%
  filter(Direction %in% c("Positive", "Non-significant")) %>%
  group_by(Direction) %>%
  select(ends_with("sum"), ends_with("sum2")) %>%
  summarize_all(list(sum)) %>%
  arrange(desc(Direction))%>%
  g(chromStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:15) %>%
  mutate(Type = "Positive")

bothC <- chromSum %>%
  filter(Direction %in% c("Both", "Non-significant")) %>%
  group_by(Direction) %>%
  select(ends_with("sum"), ends_with("sum2")) %>%
  summarize_all(list(sum)) %>%
  g(chromStates) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:15) %>%
  mutate(Type = "Both")


png("EnrichChrom.png", width = 24.24, height = 11.95, units = "cm", res = 300)
rbind(sigC, negC, posC, bothC) %>%
  spread(par, Value) %>%
  mutate(p.val.thres = ifelse(p.val > 0.05, "P > 0.05", 
                              ifelse(p.val < 1e-3, "P < 0.001", "P < 0.05")),
         Type = factor(Type, levels = c("Significant", "Inverse", "Positive", "Both")),
         Group = factor(ifelse(Region %in% c("TssA", "TssAFlnk"), "TssProxProm",
                               ifelse(Region %in% c("Tx", "TxWk"), "ActTrans", 
                                      ifelse(Region %in% c("Enh", "EnhG"), "Enhancer", 
                                             ifelse(Region %in% c("TssBiv", "BivFlnk", "EnhBiv"), "BivReg", ifelse(Region %in% c("ReprPC", "ReprPCWk"), "ReprPoly", Region))))), 
                        levels = c("TssProxProm", "TxFlnk", "ActTrans", "Enhancer", "ZNF.Rpts", "BivFlnk", "BivReg", "Het", "ReprPoly", "Quies"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "") +
  facet_grid(~ Group, scales = "free", space = "free_x") +
  theme_bw() +
  theme(plot.title = element_text(family = "Calibri", face = "bold", size = 35),
        axis.title = element_text(family = "Calibri", size = 25), 
        axis.text = element_text(family = "Calibri", size = 15),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(family = "Calibri", size = 25), 
        legend.text = element_text(family = "Calibri", size = 15), 
        strip.background = element_blank(),
        strip.text.x = element_blank())
dev.off()


# Models Comparison ####
library(dplyr)
library(ggplot2)

load("results/preprocessFiles/gexpAnnotation.Rdata")

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df

load("results/MethComBatExpResidualsCellAdj/allres_simP_cpgs.Rdata")
modC <- df

mergeTB <- modU %>%
  left_join(modC, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  filter(sigPair.x == TRUE | sigPair.y == TRUE) %>%
  mutate(sigType = ifelse(sigPair.x == TRUE, ifelse(sigPair.y == TRUE, "Both", "Crude"), "Cell"))

## P-value comparison ####
png("PvalsComp.png", width = 24.79, height = 14, units = "cm", res = 300)
ggplot(mergeTB, aes(x = -log10(p.value.y), y = -log10(p.value.x), col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Not cell adjusted") + 
  ggtitle("-log10 p-values comparative") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "")  +
  theme(plot.title = element_text(family = "Calibri", size = 30),
        axis.title = element_text(family = "Calibri", size = 25), 
        axis.text = element_text(family = "Calibri", size = 15),
        legend.title = element_text(family = "Calibri", size = 25), 
        legend.text = element_text(family = "Calibri", size = 15), 
        strip.text.x = element_text(family = "Calibri", size = 25))
dev.off()

## Estimates comparison ####
mergeTB <- mergeTB %>%
  mutate(Diff.pval = pnorm(abs(FC.x - FC.y)/sqrt(SD.x**2 + SD.y**2), lower.tail = FALSE), 
         isDiff = ifelse(Diff.pval < 0.05, "Different", "Equal"))

png("EstimatesComp.png", width = 24.79, height = 14, units = "cm", res = 300)
ggplot(mergeTB, aes(x = FC.y, y = FC.x, col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Not cell adjusted") + 
  ggtitle("Estimates comparative") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "") +
  theme(plot.title = element_text(family = "Calibri", size = 30),
        axis.title = element_text(family = "Calibri", size = 25), 
        axis.text = element_text(family = "Calibri", size = 15),
        legend.title = element_text(family = "Calibri", size = 25), 
        legend.text = element_text(family = "Calibri", size = 15), 
        strip.text.x = element_text(family = "Calibri", size = 25))
dev.off()


