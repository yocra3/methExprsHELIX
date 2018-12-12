###############################################################################
## Make figures for JBI
###############################################################################

# Descriptive of CpG-TC significant pairs ####
#'#################################################################################
load("selectedPairs.Rdata")

library(ggplot2)
library(cowplot)

## CpGs
### Total CpGs
length(unique(selPairs$CpG))
# [1] 11797

### TCs x CpG
cpgSum <- data.frame(table(table(selPairs$CpG)))
cpgSum2 <- data.frame(groups = c("1", "2", "3-5", "6-10", "+10"),
                      nums = tapply(cpgSum$Freq, rep(1:5, c(1, 1, 3, 5, 9)), sum))
cpgSum2$Prop <- cpgSum2$nums/sum(cpgSum2$nums)

cpg <- ggplot(cpgSum2, aes(x = groups, y = Prop*100)) + geom_bar(stat = "identity") +
  theme_bw() + 
  scale_x_discrete(limits  = c("1", "2", "3-5", "6-10", "+10"), name = "TCs x CpG") + 
  scale_y_continuous(name = "Prop of CpGs (%)", limits = c(0, 65)) +
  theme(axis.title = element_text(family = "Myriad", size = 30), 
        axis.text = element_text(family = "Myriad", size = 25),
        strip.text.x = element_text(family = "Myriad", size = 30),
        plot.margin = unit(c(3,3,3,3), "lines"))


## TCs
### Total TCs
length(unique(selPairs$TC))
# [1] 3743

### CpG x TC
TCSum <- data.frame(table(table(selPairs$TC)))
TCSum2 <- data.frame(groups = c("1", "2", "3-5", "6-10", "+10"),
                      nums = tapply(TCSum$Freq, rep(1:5, c(1, 1, 3, 5, 53)), sum))
TCSum2$Prop <- TCSum2$nums/sum(TCSum2$nums)

tc <- ggplot(TCSum2, aes(x = groups, y = Prop*100)) + 
  geom_bar(stat = "identity", color = "darkolivegreen", fill = "darkolivegreen") +
  theme_bw() + 
  scale_x_discrete(limits  = c("1", "2", "3-5", "6-10", "+10"), name = "CpGs x TC") + 
  scale_y_continuous(name = "Proportion of TCs (%)", limits = c(0, 65)) +
  theme(axis.title = element_text(family = "Myriad", size = 30), 
        axis.text = element_text(family = "Myriad", size = 25),
        strip.text.x = element_text(family = "Myriad", size = 30),
        plot.margin = unit(c(3,3,3,3), "lines"))


png("CpGTransDescriptives.png", width = 4946.662, height = 3506.662, res = 300)
plot_grid(cpg, tc, ncol = 2)
dev.off()

# Association CpG-TC vs distance ####
#'#################################################################################
load("selectedPairs.Rdata")

library(ggplot2)

## Distance plot ####
selPairs$Sig <- ifelse(selPairs$Estimate > 0, "Positive", "Negative")

p1 <- ggplot(selPairs, aes(x = dist, fill = Sig)) + 
  geom_density(alpha = 0.5) + 
  scale_x_continuous(limits = c(-5e5, 5e5), name = "",
                     breaks = seq(-5e5, 5e5, 1e5),
                     labels = c("-500Kb","-400Kb","-300Kb", "-200Kb", "-100Kb", "0", "100Kb", "200Kb", "300Kb", "400Kb","500Kb")) +
  scale_y_continuous(name = "") + 
  scale_fill_discrete(name = "Effect") +
  theme_bw() +
  theme(axis.title = element_text(family = "Myriad", size = 30), 
        axis.text = element_text(family = "Myriad", size = 25, vjust = 0.5, angle = 90),
        legend.key.size = unit(3, 'lines'),
        legend.title = element_text(family = "Myriad", size = 30),
        legend.text = element_text(family = "Myriad", size = 25),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        strip.text.x = element_text(family = "Myriad", size = 30))

jpeg("DistanceCpGTranscript.jpg", width = 5946.662, height = 3506.662, res = 300)
p1
dev.off()

# Enrichment CpGs ####
#'#################################################################################
load("selectedPairs.Rdata")
load("../cpganno.RData")

cpganno$sel <- ifelse(cpganno$Name %in% selPairs$CpG, "meQTL", "non-meQTL")

## Genic vs intergenic
cpganno$inter <- ifelse(cpganno$UCSC_RefGene_Name == "", "Intergenic", "Genic")

t <- table(cpganno$sel, cpganno$inter)
chisq.test(t)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  table(cpganno$sel, cpganno$inter)
# X-squared = 53.444, df = 1, p-value = 2.661e-13

# OR
t[1,1]*t[2,2]/(t[1,2]*t[2,1])
[1] 1.186881
  
## Gene positions
genTab <- lapply(colnames(cpganno)[4:9], function(name){
  t <- table(cpganno$sel, ifelse(cpganno[, name], "IN", "OUT"))
})
names(genTab) <- colnames(cpganno)[4:9]
sapply(genTab, function(x) signif(chisq.test(x)$p.value, 3))
# TSS200   TSS1500      UTR5 FirstExon      Body      UTR3
# 1.79e-08  2.61e-41  1.00e+00  3.51e-12  1.50e-03  3.73e-01

## OR
sapply(genTab, function(t) t[1,1]*t[2,2]/(t[1,2]*t[2,1]))
# TSS200   TSS1500      UTR5 FirstExon      Body      UTR3
# 0.8485928 1.3509003 1.0002576 0.7719207 1.0631959 0.9574667

cpganno$inter <- ifelse(cpganno$UCSC_RefGene_Name == "", "Intergenic", "Genic")

## Chromatin states
chrom <- read.csv("../crom15.csv")
rownames(chrom) <- chrom$HT12v4.ArrayAddress

selChrom <- c("TssA", "TssAFlnk", "Enh", "EnhG", "ReprPC", "ReprPCWk", "Quies")
cpganno <- cbind(cpganno, chrom[cpganno$Name, selChrom ])

chromTab <- lapply(selChrom, function(name){
  t <- table(cpganno$sel, ifelse(cpganno[, name], "IN", "OUT"))
})
names(chromTab) <- selChrom
sapply(chromTab, function(x) signif(chisq.test(x)$p.value, 3))
# TssA TssAFlnk      Enh     EnhG   ReprPC ReprPCWk    Quies
# 4.12e-40 0.00e+00 0.00e+00 2.02e-73 1.44e-02 1.05e-03 8.84e-04


## OR
sapply(chromTab, function(t) t[1,1]*t[2,2]/(t[1,2]*t[2,1]))
# TSS200   TSS1500      UTR5 FirstExon      Body      UTR3
# 0.8485928 1.3509003 1.0002576 0.7719207 1.0631959 0.9574667

df <- data.frame(OR = sapply(chromTab, function(t) t[1,1]*t[2,2]/(t[1,2]*t[2,1])),
                 state = selChrom,
                 type = rep(c("Promoter", "Enhancer", "Quiescent"), c(2, 2, 3)))

p2 <- ggplot(df, aes(x = state, y = OR, fill = type)) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous(trans = "log2") + 
  scale_x_discrete(limits = selChrom, name = "Chromatin states") +
  theme_bw() + scale_fill_discrete(name = "", limits = c("Promoter", "Enhancer", "Quiescent")) + 
  theme(axis.title = element_text(family = "Myriad", size = 30), 
        axis.text = element_text(family = "Myriad", size = 25, vjust = 0.5, angle = 90),
        legend.key.size = unit(3, 'lines'),
        legend.title = element_text(family = "Myriad", size = 30),
        legend.text = element_text(family = "Myriad", size = 25),
        strip.text.x = element_text(family = "Myriad", size = 30))


png("ChromatinStates.png", width = 4946.662, height = 3506.662, res = 300)
p2
dev.off()
