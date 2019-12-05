#'##############################################################################
#' eQTMs in boys and girls 
#' This script contains the code for
#'  - eQTMs in sexual chromosomes
#'  - Differences between boys and girls in autosomes
#'  - Compare Non-cell adjusted  (reference) vs cell adjusted models
#'##############################################################################

## Load libraries ####
library(ggplot2)
library(cowplot)
library(tidyr)
library(UpSetR)
library(S4Vectors)
library(dplyr)


## Load datasets ####
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/gexpAnnotation.Rdata")

## Change name when loading results
load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modA <- df
featsA <- featStatsDF

load("results/MethComBatExpResidualsNoCellAdjStrat/male/allres_simP_cpgs.Rdata")
modB <- df
featsB <- featStatsDF

load("results/MethComBatExpResidualsNoCellAdjStrat/female/allres_simP_cpgs.Rdata")
modG <- df
featsG <- featStatsDF

## Get annotation
codingTCs <- subset(expAnnot, Coding == "coding")$transcript_cluster_id

### Add annotation
modBAnn <- left_join(modB, overDF, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  mutate(chr = substring(TC, 3, 4),
         chr = gsub("^0", "", chr))

modGAnn <- left_join(modG, overDF, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  mutate(chr = substring(TC, 3, 4),
         chr = gsub("^0", "", chr))

## Merge results in one tibble
mergeAll <- full_join(modG, modB, by = c("CpG", "TC")) %>%
  full_join(modA, by = c("CpG", "TC")) %>%
  as_tibble()

## Select pairs significant in at least one dataset
### Test differences in coefficients between the models
mergeAllsig <- mergeAll %>%
  filter(sigPair | sigPair.x | sigPair.y) %>%
  mutate(FC.a_x = FC.x - FC, 
         FC.a_y = FC.y - FC, 
         FC.x_y = FC.y - FC.x, 
         pval.a_x = pnorm(abs(FC - FC.x)/sqrt(SD**2 + SD.x**2), lower.tail = FALSE),
         pval.a_y = pnorm(abs(FC - FC.y)/sqrt(SD**2 + SD.y**2), lower.tail = FALSE),
         pval.x_y = pnorm(abs(FC.x - FC.y)/sqrt(SD.x**2 + SD.y**2), lower.tail = FALSE),
         fdr.a_x = p.adjust(pval.a_x, method = "BH"),
         fdr.a_y = p.adjust(pval.a_y, method = "BH"),
         fdr.x_y = p.adjust(pval.x_y, method = "BH"),
         sig = ifelse(sigPair.x, ifelse(sigPair.y, "Common", "Girls"), 
                      ifelse(sigPair.y, "Boys", "Common")),
         sig = factor(sig, levels = c("Common", "Girls", "Boys")))

## No pairs have different coefficients between male/female and using all samples
sum(mergeAllsig$fdr.a_x < 0.05, na.rm = T)
# [1] 0
sum(mergeAllsig$fdr.a_y < 0.05, na.rm = T)
# [1] 0

#### All chr
sum(mergeAllsig$fdr.x_y < 0.05, na.rm = T)
# [1] 3761

#### Autosomes
sum(mergeAllsig$fdr.x_y < 0.05 & !is.na(mergeAllsig$sigPair), na.rm = T)
# [1] 3489

## Hits in autosome chromosomes ####
## Girls - all
### Pairs
mergeAllsig %>% filter(sigPair.x & !is.na(sigPair)) %>%
  summarize(n = n())

### TCs
mergeAllsig %>% filter(sigPair.x & !is.na(sigPair)) %>%
  select(TC) %>%
  distinct() %>%
  dplyr::summarize(n = n())
### Coding TCs
mergeAllsig %>% filter(sigPair.x & !is.na(sigPair)) %>%
  filter(TC %in% codingTCs) %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
### CpGs
mergeAllsig %>% filter(sigPair.x & !is.na(sigPair)) %>%
  select(CpG) %>%
  distinct() %>%
  summarize(n = n())

## Girls specific
### Pairs
mergeAllsig %>% filter(sigPair.x & !sigPair.y & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  summarize(n = n())
### TCs
mergeAllsig %>% filter(sigPair.x & !sigPair.y & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  select(TC) %>%
  distinct() %>%
  dplyr::summarize(n = n())
### Coding TCs
mergeAllsig %>% filter(sigPair.x & !sigPair.y & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  filter(TC %in% codingTCs) %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
### CpGs
mergeAllsig %>% filter(sigPair.x & !sigPair.y & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  select(CpG) %>%
  distinct() %>%
  summarize(n = n())


## Boys - all
### Pairs
mergeAllsig %>% filter(sigPair.y & !is.na(sigPair)) %>%
  summarize(n = n())
### TCs
mergeAllsig %>% filter(sigPair.y & !is.na(sigPair)) %>%
  select(TC) %>%
  distinct() %>%
  dplyr::summarize(n = n())
### Coding TCs
mergeAllsig %>% filter(sigPair.y & !is.na(sigPair)) %>%
  filter(TC %in% codingTCs) %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
### CpGs
mergeAllsig %>% filter(sigPair.y & !is.na(sigPair)) %>%
  select(CpG) %>%
  distinct() %>%
  summarize(n = n())

## Boys specific
### Pairs
mergeAllsig %>% filter(sigPair.y & !sigPair.x & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  summarize(n = n())
### TCs
mergeAllsig %>% filter(sigPair.y & !sigPair.x & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  select(TC) %>%
  distinct() %>%
  dplyr::summarize(n = n())
### Coding TCs
mergeAllsig %>% filter(sigPair.y & !sigPair.x & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  filter(TC %in% codingTCs) %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
### CpGs
mergeAllsig %>% filter(sigPair.y & !sigPair.x & !is.na(sigPair) & fdr.x_y < 0.05) %>%
  select(CpG) %>%
  distinct() %>%
  summarize(n = n())




## Compare Estimates ####
## Compare p-values (Only pairs significant in males or females)
png("paper/Comp_sex_P_values.png", width = 3500, height = 3500, res = 300)
mergeAllsig %>%
  filter(!is.na(p.value) & (sigPair.x | sigPair.y)) %>%
  ggplot(aes(x = -log10(p.value.x), y = -log10(p.value.y), color = sig)) +
  geom_point() +
  scale_x_continuous(name = "Girls", limits = c(0, 135)) + 
  scale_y_continuous("Boys", limits = c(0, 135)) + 
  ggtitle("-log10 p-values comparative") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "", values = c("#777777", "#8900f9", "#00c4aa")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
dev.off()


## Compare FC (Only pairs significant in males or females)
top <- mergeAllsig %>%
  filter(!is.na(p.value) & sig == "Common") %>%
  ggplot(aes(x = FC.x/10, y = FC.y/10)) +
  geom_point(color = "#777777") +
  scale_x_continuous(name = "Girls") + 
  scale_y_continuous("Boys") + 
  ggtitle("FC comparative (shared eQTMs)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = "lm", color = "#777777")

bottom <- mergeAllsig %>%
  filter(!is.na(p.value) & sig != "Common") %>%
  ggplot(aes(x = FC.x/10, y = FC.y/10, color = sig)) +
  geom_point() +
  scale_x_continuous(name = "Girls") + 
  scale_y_continuous("Boys") + 
  ggtitle("FC comparative (specific eQTMs)") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_discrete(name = "") +
  facet_wrap(. ~ sig) +
  scale_color_manual(name = "", values = c("#8900f9", "#00c4aa")) +
  geom_smooth(method = "lm") 


png("paper/Comp_sex_FC.png", width = 3500, height = 3500, res = 300)
plot_grid(top, bottom, nrow = 2)
dev.off()

filter(mergeAllsig, !is.na(p.value) & sig == "Common") %>% lm(FC.y ~ FC.x, .) %>% summary()
filter(mergeAllsig, !is.na(p.value) & sig == "Girls") %>% lm(FC.y ~ FC.x, .) %>% summary()
filter(mergeAllsig, !is.na(p.value) & sig == "Boys") %>% lm(FC.y ~ FC.x, .) %>% summary()


## Examples
filter(mergeAllsig, !is.na(p.value) & sig == "Girls") 

## Hits in sexual chromosomes ####
## Girls ####
### Pairs
modGAnn %>% filter(sigPair & chr == "X") %>%
  summarize(n = n())
# 156
### TCs
modGAnn %>% filter(sigPair & chr == "X") %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
# 71
### Coding TCs
modGAnn %>% filter(sigPair & chr == "X") %>%
  filter(TC %in% codingTCs) %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
# 49
### CpGs
modGAnn %>% filter(sigPair & chr == "X") %>%
  select(CpG) %>%
  distinct() %>%
  summarize(n = n())
# 113

## Boys ####
### Chr X
### Pairs
modBAnn %>% filter(sigPair & chr == "X") %>%
  summarize(n = n())
# 476
### TCs
modBAnn %>% filter(sigPair & chr == "X") %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
# 170
### Coding TCs
modBAnn %>% filter(sigPair & chr == "X") %>%
  filter(TC %in% codingTCs) %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
# 130
### CpGs
modBAnn %>% filter(sigPair & chr == "X") %>%
  select(CpG) %>%
  distinct() %>%
  summarize(n = n())
# 327

### Chr Y
### Pairs
modBAnn %>% filter(sigPair & chr == "Y") %>%
  summarize(n = n())
# 23
### TCs
modBAnn %>% filter(sigPair & chr == "Y") %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
# 8
### Coding TCs
modBAnn %>% filter(sigPair & chr == "Y") %>%
  filter(TC %in% codingTCs) %>%
  select(TC) %>%
  distinct() %>%
  summarize(n = n())
# 2
### CpGs
modBAnn %>% filter(sigPair & chr == "Y") %>%
  select(CpG) %>%
  distinct() %>%
  summarize(n = n())
# 15

## Distribution CpGs/TC ####
modGAnn %>% 
  filter(sigPair & chr == "X") %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   1.000   2.197   2.000  18.000
modBAnn %>% 
  filter(sigPair & chr == "X") %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.0     1.0     1.0     2.8     3.0    29.0

modBAnn %>% 
  filter(sigPair & chr == "Y") %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.750   2.500   2.875   4.250   5.000


modGAnn %>% 
  filter(sigPair & chr == "X") %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   1.000   1.381   1.000   7.000
modBAnn %>% 
  filter(sigPair & chr == "X") %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   1.000   1.456   2.000   8.000

modBAnn %>% 
  filter(sigPair & chr == "Y") %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   2.000   1.533   2.000   2.000



png("paper/eQTMs_sexSpecific_Chr_distr.png", width = 2000, height = 1500, res = 300)
mergeAllsig %>% 
  filter(p.adjust(mergeAllsig$pval.x_y, "BH") < 0.05) %>%
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


mergeSexSp <- mergeAllsig %>% 
  filter(p.adjust(mergeAllsig$pval.x_y, "BH") < 0.05) 
  
mergeSexSp %>%
  summarize(Girls = sum(sigPair.x & !sigPair.y),
            GirlsAll = sum(sigPair.x & !sigPair, na.rm = TRUE),
            Boys = sum(sigPair.y & !sigPair.x),
            BoysAll = sum(sigPair.y & !sigPair, na.rm = TRUE))
# Girls GirlsAll  Boys BoysAll
# <int>    <int> <int>   <int>
#  1159      924  2432    1294

## Girls
mergeSexSp %>%
  filter(sigPair.x & !sigPair.y) %>%
  summarize(TCs = length(unique(TC)),
            TCcoding = sum(unique(TC) %in% codingTCs),
            CpGs = length(unique(CpG)))
# TCs TCcoding  CpGs
# <int>    <int> <int>
#  932      565  1093

## Girls sp
mergeSexSp %>%
  filter(sigPair.x & !sigPair) %>%
  summarize(TCs = length(unique(TC)),
            TCcoding = sum(unique(TC) %in% codingTCs),
            CpGs = length(unique(CpG)))
# TCs TCcoding  CpGs
# <int>    <int> <int>
#   778      447   885

## Boys
mergeSexSp %>%
  filter(sigPair.y & !sigPair.x) %>%
  summarize(TCs = length(unique(TC)),
            TCcoding = sum(unique(TC) %in% codingTCs),
            CpGs = length(unique(CpG)))
# TCs TCcoding  CpGs
# <int>    <int> <int>
#  1484     1062  2117

## Boys sp
mergeSexSp %>%
  filter(sigPair.y & !sigPair) %>%
  summarize(TCs = length(unique(TC)),
            TCcoding = sum(unique(TC) %in% codingTCs),
            CpGs = length(unique(CpG)))
# TCs TCcoding  CpGs
# <int>    <int> <int>
#  1020      677  1220

mergeSexSp  %>% 
  filter(sigPair.y & sigPair.x) %>% 
  mutate(cot = sign(FC.x) == sign(FC.y)) %>% 
  filter(!cot) %>% 
  select(CpG, TC, FC.x, p.value.x, FC.y, p.value.y)
# A tibble: 2 x 6
# CpG        TC               FC.x   p.value.x   FC.y p.value.y
# <chr>      <chr>           <dbl>       <dbl>  <dbl>     <dbl>
# cg12026625 TC0X000267.hg.1  2.24 0.000000506 -94.1   8.73e-33
# cg27088126 TC0X001945.hg.1  1.59 0.00000504   -2.82  6.85e- 5


## Gene Enrichment ####
## Copy Functions from eQTM_interpretation.R
runGO <- function(genes, df){
  selGenes <- df %>%
    as_tibble() %>%
    dplyr::select(TC) %>%
    distinct() %>%
    mutate(sig = factor(ifelse(TC %in% genes, 1, 0)))
  go <- computeGOs(selGenes)
  go
}


## Subtypes
goSex <- list(SexDifs = runGO(unique(mergeSexSp$TC), modG),
              GirlSp = runGO(unique(filter(mergeSexSp, sigPair.x & !sigPair.y)$TC), modG),
              BoySp = runGO(unique(filter(mergeSexSp, sigPair.y & !sigPair.x)$TC), modG)
)

save(goSex, file = "paper/sexSpGOobjects.Rdata")

## It does not work on server. Run locally.
library(GOfuncR)
library(dplyr)
server <- "//isg10174/data/WS_HELIX/HELIX_analyses/expr_met_SM/paper/"

load(paste0(server, "sexSpGOobjects.Rdata"))

addImmunityInfo <- function(tab){
  get_parent_nodes(tab$GO.ID) %>%
    group_by(child_go_id) %>%
    summarize(tag = ifelse(any(c("adaptive immune response", "lymphocyte activation", "antigen processing and presentation") %in% parent_name), "adaptive",
                           ifelse(any(c("innate immune response", "myeloid leukocyte activation", "myeloid leukocyte mediated immunity", "myeloid leukocyte differentiation", "myeloid leukocyte cytokine production") %in% parent_name), "innate",
                                  ifelse(any(c("response to cytokine", "immune system process") %in% parent_name), "general", "none"))),
              go = parent_name[1]) %>%
    mutate(immune = tag != "none", 
           GO.ID = child_go_id) %>%
    right_join(tab)  %>%
    select(GO.ID, go, w0, classic, immune, tag)
}

sexTabs <- lapply(goSex, function(x){
  addImmunityInfo(x$tab)
})


lapply(names(sexTabs), function(x){
  write.table(sexTabs[[x]][, c("GO.ID", "go", "w0", "classic", "tag")], 
              file = paste0(server, "/GOs", x, ".txt"), 
              quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
  
})
lapply(sexTabs, function(x){
  c(n = nrow(x), imm = sum(x$immune), immP = round(mean(x$immune)*100, 1), 
    adap = sum(x$tag == "adaptive"), innate = sum(x$tag == "innate"))
})
