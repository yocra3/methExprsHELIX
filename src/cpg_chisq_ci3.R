
## Setup

setwd("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM")
load("res2.data.RData")
head(res2.data)
bonf <- 0.05/nrow(res2.data)
bonf

res2.data$bonf <- ifelse(res2.data$p.value < bonf, "sig", "nosig")
rsig <- res2.data[res2.data$bonf == "sig", ]

rsig <- rsig[!duplicated(rsig$CpG), ]
dim(rsig) # 8907   21

rnosig <- res2.data[!(res2.data$CpG %in% rsig$CpG), ]
rnosig <- rnosig[!(duplicated(rnosig$CpG)), ]
dim(rnosig) # 377511     21

r <- rbind(rsig, rnosig)
dim(r) # 386418     21


## Genomic position

r$g <- ifelse(r$UCSC_RefGene_Name != "", "Genic", "Intergenic")

g <- table(r$bonf, r$g)
g <- as.data.frame.matrix(g)
g <- g[c(2,1),]
g

chisq.test(g) # 1.818e-14
or <- (g[1, 1] * g[2, 2]) / (g[1, 2] * g[2, 1])
cat("Odds Ratio: ", or, "\n", sep = "") # 1.231447
print(exp(log(or)+1.96*sqrt(sum(1/g)))) # 1.298803
print(exp(log(or)-1.96*sqrt(sum(1/g)))) # 1.167584


## Relation to island

rti <- table(r$bonf, r$Relation_to_Island)
rti <- as.data.frame.matrix(rti)
rti <- rti[c(2,1),]
rti

coln <- NULL
pvals1 <- NULL
h1 <- NULL
l1 <- NULL
or1 <- NULL
for (i in 1:6) {
  coln <- c(coln, colnames(rti)[i])
  df <- cbind(rti[,i],rowSums(rti[,-i]))
  xt <- chisq.test(df)
  or <- (df[1, 1] * df[2, 2]) / (df[1, 2] * df[2, 1])
  or1 <- c(or1, or)
  h1 <- c(h1, exp(log(or)+1.96*sqrt(sum(1/df))))
  l1 <- c(l1, exp(log(or)-1.96*sqrt(sum(1/df))))
  pvals1 <- c(pvals1, xt$p.value)
}

island <- data.frame("Island" = coln,
                    "Sig_vs_nosig_pval" = pvals1,
                    "OR" = or1,
                    "CIH" = h1,
                    "CIL" = l1)
island

#    Island Sig_vs_nosig_pval        OR       CIH       CIL
# 1  Island      4.523043e-61 0.6631389 0.6965287 0.6313498
# 2 N_Shelf      1.099837e-04 0.8108555 0.9013254 0.7294665
# 3 N_Shore      1.266199e-62 1.5715964 1.6578614 1.4898200
# 4 OpenSea      1.805701e-06 0.8962596 0.9374013 0.8569236
# 5 S_Shelf      6.874332e-03 0.8590976 0.9582091 0.7702376
# 6 S_Shore      1.657955e-77 1.7110011 1.8115236 1.6160566


## Relative position

coln <- NULL
pvals1 <- NULL
h1 <- NULL
l1 <- NULL
or1 <- NULL
for (i in 14:19) {
  cat("\n########################################################\n\n")
  print(colnames(r)[i])
  rp <- table(r$bonf, r[,i])
  rp <- as.data.frame.matrix(rp)
  rp <- rp[,c(2,1)]
  coln <- c(coln, colnames(r)[i])
  
  cat("\n## Sig vs nosig ##\n\n")
  rp1 <- rbind(nosig = rp[1, ], sig = colSums(rp[-1, ]))
  rp1 <- rp1[c(2,1),]
  print(rp1)
  xt <- chisq.test(rp1)
  print(xt)
  or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])
  cat("Odds Ratio: ", or, "\n", sep = "")
  or1 <- c(or1, or)
  h1 <- c(h1, exp(log(or)+1.96*sqrt(sum(1/rp1))))
  l1 <- c(l1, exp(log(or)-1.96*sqrt(sum(1/rp1))))
  pvals1 <- c(pvals1, xt$p.value)

}

relp <- data.frame("RelPosition" = coln,
                    "Sig_vs_nosig_pval" = pvals1,
                    "OR" = or1,
                    "CIH" = h1,
                    "CIL" = l1)
print(relp)

#   RelPosition Sig_vs_nosig_pval        OR       CIH       CIL
# 1      TSS200      3.804229e-04 0.8894715 0.9485986 0.8340298
# 2     TSS1500      6.865677e-36 1.3748566 1.4453984 1.3077575
# 3        UTR5      1.241552e-01 1.0482036 1.1123781 0.9877315
# 4   FirstExon      6.901174e-08 0.7961830 0.8648321 0.7329832
# 5        Body      1.001492e-01 1.0373011 1.0832756 0.9932777
# 6        UTR3      9.485876e-02 0.9096585 1.0148187 0.8153954


## Chromatin States

library("minfi")
load("gset_sm_1.RData") # gset
cs <- rowData(gset)
cs <- as.data.frame(cs[, c(1,29:43)])
colnames(cs)[1] <- "CpG"
dim(cs) # 386518     16
head(cs)
r <- merge(r, cs, by = "CpG", sort = FALSE)

coln <- NULL
pvals1 <- NULL
h1 <- NULL
l1 <- NULL
or1 <- NULL
for (i in 23:37) {
  cat("\n########################################################\n\n")
  print(colnames(r)[i])
  rp <- table(r$bonf, r[,i])
  rp <- as.data.frame.matrix(rp)
  rp <- rp[,c(2,1)]
  coln <- c(coln, colnames(r)[i])
  
  cat("\n## Sig vs nosig ##\n\n")
  rp1 <- rbind(nosig = rp[1, ], sig = colSums(rp[-1, ]))
  rp1 <- rp1[c(2,1),]
  print(rp1)
  xt <- chisq.test(rp1)
  print(xt)
  or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])
  cat("Odds Ratio: ", or, "\n", sep = "")
  or1 <- c(or1, or)
  h1 <- c(h1, exp(log(or)+1.96*sqrt(sum(1/rp1))))
  l1 <- c(l1, exp(log(or)-1.96*sqrt(sum(1/rp1))))
  pvals1 <- c(pvals1, xt$p.value)
}

chroms <- data.frame("Chromatin_State" = coln,
                    "Sig_vs_nosig_pval" = pvals1,
                    "OR" = or1,
                    "CIH" = h1,
                    "CIL" = l1)

print(chroms)

#    Chromatin_State Sig_vs_nosig_pval        OR       CIH       CIL
# 1             TssA      7.116479e-57 1.4102799 1.4716017 1.3515133
# 2         TssAFlnk      0.000000e+00 2.3560387 2.4600495 2.2564255
# 3           TxFlnk      9.075034e-01 0.9960869 1.0562110 0.9393853
# 4             TxWk      5.751097e-60 1.4202810 1.4816126 1.3614882
# 5               Tx      9.902663e-04 0.9045750 0.9599800 0.8523677
# 6             EnhG      7.249001e-58 1.6994879 1.8144931 1.5917720
# 7              Enh      0.000000e+00 2.9000354 3.0261077 2.7792155
# 8         ZNF.Rpts      2.327439e-68 2.4674383 2.7392866 2.2225684
# 9              Het      5.115275e-02 1.0817202 1.1696392 1.0004098
# 10          TssBiv      5.237444e-05 1.1486417 1.2281237 1.0743035
# 11         BivFlnk      5.237444e-05 1.1486417 1.2281237 1.0743035
# 12          EnhBiv      4.497477e-43 1.4158204 1.4879931 1.3471484
# 13          ReprPC      8.527877e-01 1.0046641 1.0521804 0.9592936
# 14        ReprPCWk      3.665152e-08 0.8882686 0.9264772 0.8516357
# 15           Quies      2.998602e-02 1.0486771 1.0944414 1.0048265

### PLOTS ###
library(ggplot2)
t_shift <- scales::trans_new("shift",
                             transform = function(x) {x-1},
                             inverse = function(x) {x+1})

# Island

v <- c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "OpenSea")
island$Island <- factor(as.character(island$Island), levels = v)


pdf("1_enrichment_island_2.pdf", height = 6, width = 5)
ggplot(island, aes(x = Island, y = OR)) + 
  geom_bar(stat="identity", fill = "steelblue1", width = 0.5) +
  geom_errorbar(aes(ymin=CIL, ymax=CIH), width=.15) +
  scale_y_log10(limits = c(0.5, 2), breaks = c(0.1, 0.5, 1, 1.5, 2)) +
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Relation to island")
dev.off()

# Relative position

relp$RelPosition <- c("TSS200","TSS1500","5'UTR","1stExon","Body","3'UTR")
relp <- relp[c(2,1,3,4,5,6),]
relp$RelPosition <- factor(relp$RelPosition, levels = (as.character(relp$RelPosition)))

pdf("1_enrichment_relp_2.pdf", height = 6, width = 5)
ggplot(relp, aes(x = RelPosition, y = OR)) + 
  geom_bar(stat="identity", fill = "steelblue1", width = 0.5) +
  geom_errorbar(aes(ymin=CIL, ymax=CIH), width=.15) +
  scale_y_log10(limits = c(0.5, 2), breaks = c(0.1, 0.5, 1, 1.5, 2)) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1) +
  xlab("Relative position")
dev.off()


# Chromatin States
chroms <- chroms[c(1:3,5,4,6:15),]
x <- chroms$Chromatin_State
chroms$Chromatin_State <- factor(x, levels = as.character(x))

pdf("1_enrichment_cs.pdf", height = 6, width = 8)
ggplot(chroms, aes(x = Chromatin_State, y = OR)) + 
  geom_bar(stat="identity", fill = "steelblue1", width = 0.7) +
  geom_errorbar(aes(ymin=CIL, ymax=CIH), width=.15) +
  scale_y_log10(limits = c(0.1, 4), breaks = c(0.1, 0.5, 1, 2, 4)) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1) +
  xlab("Chromatin States") +
  ylab("Odds Ratio (Significance)")
dev.off()


