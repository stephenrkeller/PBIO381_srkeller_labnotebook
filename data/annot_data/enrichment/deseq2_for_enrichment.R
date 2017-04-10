setwd("~/Dropbox/1_Research/SSW/DGE")
library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

### Set up group model


colData$group <- factor(paste0(colData$location, ".", colData$health, ".", colData$score))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
dim(dds)

dds <- dds[ rowSums(counts(dds)) > 100, ]

dds <- dds[sample(nrow(dds), 1200), ]
dim(dds)

dds <- DESeq(dds, parallel=T)

resultsNames(dds)
# [1] "Intercept"    "groupint.H.0" "groupint.S.1" "groupint.S.2" "groupint.S.3" "groupint.S.4"
# [7] "groupint.S.5" "groupsub.H.0" "groupsub.S.1" "groupsub.S.2" "groupsub.S.3"

# Just focus on the intertidal animals, contrasting between healthy and early stages of sickness (score 1 or 2)
res <- results(dds, contrast=list( c("groupint.H.0"), c("groupint.S.1","groupint.S.2")), listValues=c(1/2, -1/2))
summary(res)
# 
# out of 12952 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 99, 0.76% 
# LFC < 0 (down)   : 9, 0.069% 
# outliers [1]     : 363, 2.8% 
# low counts [2]   : 1995, 15% 
# (mean count < 5)

res <- results(dds, contrast=list( c("groupint.H.0","groupsub.H.0"), c("groupint.S.1","groupsub.S.1")), listValues=c(1/2, -1/2))
summary(res)

# out of 12952 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 609, 4.7% 
# LFC < 0 (down)   : 42, 0.32% 
# outliers [1]     : 363, 2.8% 
# low counts [2]   : 1514, 12% 
# (mean count < 5)

res <- results(dds, contrast=list( c("groupint.H.0"), c("groupint.S.1")), listValues=c(1/2, -1/2))
summary(res)
# out of 12952 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 659, 5.1% 
# LFC < 0 (down)   : 36, 0.28% 
# outliers [1]     : 363, 2.8% 
# low counts [2]   : 1514, 12% 
# (mean count < 5)


res <- res[order(res$padj),]
head(res)
# log2 fold change (MAP): 0.5 groupint.H.0 vs 0.5 groupint.S.1+groupint.S.2 
# Wald test p-value: 0.5 groupint.H.0 vs 0.5 groupint.S.1+groupint.S.2 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN35473_c0_g1_TRINITY_DN35473_c0_g1_i1_g.5722_m.5722     36.83831       3.680482 0.5232499  7.033890
# TRINITY_DN40577_c0_g1_TRINITY_DN40577_c0_g1_i1_g.10012_m.10012   39.36954       3.571670 0.5277764  6.767393
# TRINITY_DN40310_c0_g1_TRINITY_DN40310_c0_g1_i1_g.9594_m.9594     53.44457       3.274379 0.4881744  6.707397
# TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710   7429.45917       6.187736 0.9402100  6.581228
# TRINITY_DN35868_c0_g1_TRINITY_DN35868_c0_g1_i1_g.5949_m.5949     20.48137       3.251657 0.6109193  5.322564
# TRINITY_DN38692_c0_g1_TRINITY_DN38692_c0_g1_i1_g.7973_m.7973     31.90442       3.035705 0.6008663  5.052213
# pvalue             padj
# <numeric>        <numeric>
#   TRINITY_DN35473_c0_g1_TRINITY_DN35473_c0_g1_i1_g.5722_m.5722   0.000000000002008539 0.00000002128248
# TRINITY_DN40577_c0_g1_TRINITY_DN40577_c0_g1_i1_g.10012_m.10012 0.000000000013112380 0.00000006946939
# TRINITY_DN40310_c0_g1_TRINITY_DN40310_c0_g1_i1_g.9594_m.9594   0.000000000019812634 0.00000006997822
# TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710   0.000000000046657949 0.00000012359691
# TRINITY_DN35868_c0_g1_TRINITY_DN35868_c0_g1_i1_g.5949_m.5949   0.000000102314739823 0.00021682539663
# TRINITY_DN38692_c0_g1_TRINITY_DN38692_c0_g1_i1_g.7973_m.7973   0.000000436719528985 0.00077124668819




write.csv(res, file = "results_int_H0vsS1.csv", row.names = T, quote = F)

