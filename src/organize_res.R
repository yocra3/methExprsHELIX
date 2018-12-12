library("minfi")
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_preproc/gene_expression/
     annotation/HTAv2_0.na36.hg19.transcript_FINAL4.Rdata")
tcanno <- HTAv2_0.na36.hg19.transcript_FINAL4
rm(HTAv2_0.na36.hg19.transcript_FINAL4)
# load("/home/isglobal.lan/smari/data/smari/gset_sm_1.RData")
# cpganno <- rowData(gset)
# cpganno <- as.data.frame(cpganno[, c(5,9,12,23:28)])
# save(cpganno, file = "cpganno.RData")
# rm(gset)
load("/home/isglobal.lan/smari/data/smari/cpganno.RData")
load("/home/isglobal.lan/smari/data/smari/overlaps_d.RData")
load("/home/isglobal.lan/smari/data/smari/fil.RData")


### RESULTS 1
setwd("/home/isglobal.lan/smari/data/smari")
load("/home/isglobal.lan/smari/data/smari/data_res1/all_output.RData")
res1 <- all_output
rm(all_output)
res1 <- res1[!(res1$TC %in% fil), ]
res1$FDR <- p.adjust(res1$p.value, method = "fdr")
res1 <- merge(res1, tcanno[,c(1,24)], by.x = "TC",
              by.y = "transcript_cluster_id", sort = F)
save(res1, file = "res1.RData")

res1.data <- merge(res1, overlaps_d[, c(1,2,6,8)], by.x = c("CpG", "TC"),
                   by.y = c("CpG", "TC"), sort = F)
res1.data <- merge(res1.data, cpganno, by.x = "CpG", by.y = "Name", sort = F)
res1.data <- merge(res1.data, tcanno[,c(1,15)], by.x = "TC",
                                     by.y = "transcript_cluster_id", sort = F)
save(res1.data, file = "res1.data.RData")


### RESULTS 2
setwd("/home/isglobal.lan/smari/data/smari")
load("/home/isglobal.lan/smari/data/smari/data_res2/all_output.RData")
res2 <- all_output
rm(all_output)
res2 <- res2[!(res2$TC %in% fil), ]
res2$FDR <- p.adjust(res2$p.value, method = "fdr")
res2 <- merge(res2, tcanno[,c(1,24)], by.x = "TC",
              by.y = "transcript_cluster_id", sort = F)
save(res2, file = "res2.RData")

res2.data <- merge(res2, overlaps_d[, c(1,2,6,8)], by.x = c("CpG", "TC"),
                   by.y = c("CpG", "TC"), sort = F)
res2.data <- merge(res2.data, cpganno, by.x = "CpG", by.y = "Name", sort = F)
res2.data <- merge(res2.data, tcanno[,c(1,15)], by.x = "TC",
                   by.y = "transcript_cluster_id", sort = F)
save(res2.data, file = "res2.data.RData")

### RESULTS 3
setwd("/home/isglobal.lan/smari/data/smari")
load("/home/isglobal.lan/smari/data/smari/data_res3/all_output.RData")
res3 <- all_output
rm(all_output)
res3 <- res3[!(res3$TC %in% fil), ]
res3$FDR <- p.adjust(res3$p.value, method = "fdr")
res3 <- merge(res3, tcanno[,c(1,24)], by.x = "TC",
              by.y = "transcript_cluster_id", sort = F)
save(res3, file = "res3.RData")

res3.data <- merge(res3, overlaps_d[, c(1,2,6,8)], by.x = c("CpG", "TC"),
                   by.y = c("CpG", "TC"), sort = F)
res3.data <- merge(res3.data, cpganno, by.x = "CpG", by.y = "Name", sort = F)
res3.data <- merge(res3.data, tcanno[,c(1,15)], by.x = "TC",
                   by.y = "transcript_cluster_id", sort = F)
save(res3.data, file = "res3.data.RData")


load("res1.RData")
load("res2.RData")
load("res3.RData")
load("res1.data.RData")
load("res2.data.RData")
load("res3.data.RData")

### PLOTS
library(ggplot2)
setwd("/home/isglobal.lan/smari/data/smari/plots")

# width - p-value
png("pv1.png", type = "cairo", width = 1280, height = 720)
ggplot(res1.data, aes(x = dif, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(-5e5, -2.5e5, 0, 2.5e5, 5e5))
dev.off()

png("pv2.png", type = "cairo", width = 1280, height = 720)
ggplot(res2.data, aes(x = dif, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(-5e5, -2.5e5, 0, 2.5e5, 5e5))
dev.off()

png("pv3.png", type = "cairo", width = 1280, height = 720)
ggplot(res3.data, aes(x = dif, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(-5e5, -2.5e5, 0, 2.5e5, 5e5))
dev.off()

# width - estimate
png("est1.png", type="cairo", width = 1280, height = 720)
ggplot(res1.data, aes(x = abs(dif), y = Estimate)) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(0, 2.5e5, 5e5))
dev.off()

png("est2.png", type="cairo", width = 1280, height = 720)
ggplot(res2.data, aes(x = abs(dif), y = Estimate)) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(0, 2.5e5, 5e5))
dev.off()

png("est3.png", type="cairo", width = 1280, height = 720)
ggplot(res3.data, aes(x = abs(dif), y = Estimate)) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(0, 2.5e5, 5e5))
dev.off()

# estimate - p-value
png("est.pv.1.png", type = "cairo", width = 1280, height = 720)
ggplot(data = res1, aes(x = Estimate, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(0.0001754), color = "steelblue")
dev.off()

png("est.pv.2.png", type = "cairo", width = 1280, height = 720)
ggplot(data = res2, aes(x = Estimate, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(0.0001678), color = "steelblue")
dev.off()

png("est.pv.3.png", type = "cairo", width = 1280, height = 720)
ggplot(data = res3, aes(x = Estimate, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(0.0003025), color = "steelblue")
dev.off()

# width - significant p.values (FDR < 0.05) (density) (separating estimates)
dens1 <- res1.data[res1.data$FDR < 0.05, ]
dens1$effect <- as.factor(ifelse(dens1$Estimate >= 0, "Positive", "Negative"))
dens2 <- res2.data[res2.data$FDR < 0.05, ]
dens2$effect <- as.factor(ifelse(dens2$Estimate >= 0, "Positive", "Negative"))
dens3 <- res3.data[res3.data$FDR < 0.05, ]
dens3$effect <- as.factor(ifelse(dens3$Estimate >= 0, "Positive", "Negative"))

png("pvdens1.png", type = "cairo", width = 1280, height = 720)
ggplot(dens1, aes(x = abs(dif), color = effect)) +
  geom_density() +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(0, 2.5e5, 5e5))
dev.off()

png("pvdens2.png", type = "cairo", width = 1280, height = 720)
ggplot(dens2, aes(x = abs(dif),  color = effect)) +
  geom_density() +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(0, 2.5e5, 5e5))
dev.off()

png("pvdens3.png", type = "cairo", width = 1280, height = 720)
ggplot(dens3, aes(x = abs(dif),  color = effect)) +
  geom_density() +
  theme_classic(base_size = 22) +
  scale_x_continuous(breaks=c(0, 2.5e5, 5e5))
dev.off()

# sign p.values/chromosome
chr <- sapply(1:22, function(x) paste("chr", x, sep = ""))
c1 <- as.data.frame(table(res1.data[res1.data$FDR < 0.05, "chrCpG"]))
c1$Var1 <- factor(c1$Var1, levels = chr)
c1$Freq2 <- round(c1$Freq/sum(c1$Freq)*100, digits = 2)
c2 <- as.data.frame(table(res2.data[res2.data$FDR < 0.05, "chrCpG"]))
c2$Var1 <- factor(c2$Var1, levels = chr)
c2$Freq2 <- round(c2$Freq/sum(c2$Freq)*100, digits = 2)
c3 <- as.data.frame(table(res3.data[res3.data$FDR < 0.05, "chrCpG"]))
c3$Var1 <- factor(c3$Var1, levels = chr)
c3$Freq2 <- round(c3$Freq/sum(c3$Freq)*100, digits = 2)

png("chr1.png", type = "cairo", width = 1280, height = 720)
ggplot(c1, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Freq2), vjust = -1, color = "black", size = 5) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png("chr2.png", type = "cairo", width = 1280, height = 720)
ggplot(c2, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Freq2), vjust = -1, color = "black", size = 5) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png("chr3.png", type = "cairo", width = 1280, height = 720)
ggplot(c3, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Freq2), vjust = -1, color = "black", size = 5) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# p-values between results
d <- data.frame("d1" = res1$p.value,
                "d2" = res2$p.value,
                "d3" = res3$p.value)

png("pv_1_2.png", type = "cairo", height = 1080, width = 1080)
ggplot(d, aes(x = d1, y = d2)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  theme_classic(base_size = 30)
dev.off()

png("pv_1_3.png", type = "cairo", height = 1080, width = 1080)
ggplot(d, aes(x = d1, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  theme_classic(base_size = 30)
dev.off()

png("pv_2_3.png", type = "cairo", height = 1080, width = 1080)
ggplot(d, aes(x = d2, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  theme_classic(base_size = 30)
dev.off()

# estimates between results
d <- data.frame("d1" = res1$Estimate,
                "d2" = res2$Estimate,
                "d3" = res3$Estimate)

png("est_1_2.png", type = "cairo", height = 1080, width = 1080)
ggplot(d, aes(x = d1, y = d2)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  theme_classic(base_size = 30)
dev.off()

png("est_1_3.png", type = "cairo", height = 1080, width = 1080)
ggplot(d, aes(x = d1, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  theme_classic(base_size = 30)
dev.off()

png("est_2_3.png", type = "cairo", height = 1080, width = 1080)
ggplot(d, aes(x = d2, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  theme_classic(base_size = 30)
dev.off()

### Correlations
cor(res1.data$p.value, res2.data$p.value) # 0.9794996
cor(res1.data$p.value, res3.data$p.value) # 0.7427808
cor(res2.data$p.value, res3.data$p.value) # 0.7564384

cor(res1.data$Estimate, res2.data$Estimate) # 0.9934421
cor(res1.data$Estimate, res3.data$Estimate) # 0.9301162
cor(res2.data$Estimate, res3.data$Estimate) # 0.9379744

### CpG vs TC

library(minfi)
library(reshape2)

# Res 1
load("/home/isglobal.lan/smari/data/smari/summ_exp_sm_res1.RData") # se_res
load("/home/isglobal.lan/smari/data/smari/gset_sm_res1.RData") # gset_res
ex <- assays(se_res)$expr
met <- assays(gset_res)$Beta
table(colnames(ex) == colnames(met))

sig <- res1[res1$FDR < 0.05, ]
sig <- sig[sample(nrow(sig), 8), ]
ex.sig <- ex[sig$TC, ]
met.sig <- met[sig$CpG, ]
ex.sig2 <- melt(ex.sig)
met.sig2 <- melt(met.sig)
table(ex.sig2$Var2 == met.sig2$Var2)
table(sig$TC == levels(ex.sig2$Var1))
table(sig$CpG == levels(met.sig2$Var1))
d <- cbind(met.sig2, ex.sig2)
colnames(d) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")

png("sig_exp_met_res1.png", type = "cairo", width = 1080, height = 1080)
ggplot(d, aes(x = met, y = exp, color = CpG)) +
  geom_point(alpha = 0.5, size = 3) +
  theme_classic(base_size = 30)
dev.off()

nosig <- res1[res1$FDR > 0.05, ]
nosig <- nosig[sample(nrow(nosig), 8), ]
ex.nosig <- ex[nosig$TC, ]
met.nosig <- met[nosig$CpG, ]
ex.nosig2 <- melt(ex.nosig)
met.nosig2 <- melt(met.nosig)
table(ex.nosig2$Var2 == met.nosig2$Var2)
table(nosig$TC == levels(ex.nosig2$Var1))
table(nosig$CpG == levels(met.nosig2$Var1))
d2 <- cbind(met.nosig2, ex.nosig2)
colnames(d2) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")

png("nosig_exp_met_res1.png", type = "cairo", width = 1080, height = 1080)
ggplot(d2, aes(x = met, y = exp, color = CpG)) +
  geom_point(alpha = 0.5, size = 3) +
  theme_classic(base_size = 30)
dev.off()


# Res 2
load("/home/isglobal.lan/smari/data/smari/summ_exp_sm_res2.RData") # se_res
load("/home/isglobal.lan/smari/data/smari/gset_sm_1.RData") # gset
ex <- assays(se_res)$expr
met <- assays(gset)$Beta
table(colnames(ex) == colnames(met))

sig <- res2[res2$FDR < 0.05, ]
sig <- sig[sample(nrow(sig), 8), ]
ex.sig <- ex[sig$TC, ]
met.sig <- met[sig$CpG, ]
ex.sig2 <- melt(ex.sig)
met.sig2 <- melt(met.sig)
table(ex.sig2$Var2 == met.sig2$Var2)
table(sig$TC == levels(ex.sig2$Var1))
table(sig$CpG == levels(met.sig2$Var1))
d <- cbind(met.sig2, ex.sig2)
colnames(d) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")

png("sig_exp_met_res2.png", type = "cairo", width = 1080, height = 1080)
ggplot(d, aes(x = met, y = exp, color = CpG)) +
  geom_point(alpha = 0.5, size = 3) +
  theme_classic(base_size = 30)
dev.off()

nosig <- res2[res2$FDR > 0.05, ]
nosig <- nosig[sample(nrow(nosig), 8), ]
ex.nosig <- ex[nosig$TC, ]
met.nosig <- met[nosig$CpG, ]
ex.nosig2 <- melt(ex.nosig)
met.nosig2 <- melt(met.nosig)
table(ex.nosig2$Var2 == met.nosig2$Var2)
table(nosig$TC == levels(ex.nosig2$Var1))
table(nosig$CpG == levels(met.nosig2$Var1))
d2 <- cbind(met.nosig2, ex.nosig2)
colnames(d2) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")

png("nosig_exp_met_res2.png", type = "cairo", width = 1080, height = 1080)
ggplot(d2, aes(x = met, y = exp, color = CpG)) +
  geom_point(alpha = 0.5, size = 3) +
  theme_classic(base_size = 30)
dev.off()
