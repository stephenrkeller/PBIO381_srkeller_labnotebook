# Sample script provided by M. Pespeni, Feb 27, 2017
# Modified by S. Keller

library("DESeq2")

library("ggplot2")

setwd("/Users/srkeller/github/PBIO381_srkeller_labnotebook/data/DGE_data/round2")
countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

countDataint <- countData[,1:48]
countDatasub <- countData[,49:77]

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

colDataint <- subset(colData, location=="int")
colDatasub <- subset(colData, location=="sub")


#################### MODEL NUMBER 1a: TEST EFFECT OF HEALTH CONTROLLING FOR LOCATION

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)

# In this typical model, the "health effect" represents the overall effect of health status 
# controlling for differences due to location page 26-27 manual DESeq2.pdf
# The last term in the model is what is tested.  In this case health. You need to rearrange the order 
# of the factors in the model to test for the overall effect of a diff. factor.
# This is not the same as an interaction.


dim(dds)
#[1] 13053    77

dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
#[1] 12954    77  # at > 100; we loose many fewer genes


colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) #sets that "healthy is the reference

dds <- DESeq(dds) 

res <- results(dds)
res <- res[order(res$padj),]
head(res)

summary(res)

# out of 1199 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 27, 2.3% 
# LFC < 0 (down)   : 14, 1.2% 
# outliers [1]     : 31, 2.6% 
# low counts [2]   : 743, 62% 
# (mean count < 25)

Try extracting just the genes showing significant DE at a given threshold of 1% (=0.01)
res_sig <- res[which(res$padj<0.01),]
dim(res_sig)


#################### MODEL NUMBER 1int: TEST EFFECT OF HEALTH FOR JUST LOCATION=INT

dds_int <- DESeqDataSetFromMatrix(countData = countDataint, colData = colDataint, design = ~ health)

dim(dds_int)
#[1] 13053    48

dds_int <- dds_int[ rowSums(counts(dds_int)) > 100, ]
dim(dds_int)
#[1] 12405    48  # at > 100

colData(dds_int)$health <- factor(colData(dds_int)$health, levels=c("H","S")) #sets that "healthy is the reference

dds_int <- DESeq(dds_int) 

res_int <- results(dds_int)
res_int <- res_int[order(res_int$padj),]
head(res_int)

summary(res_int)

#################### MODEL NUMBER 1sub: TEST EFFECT OF HEALTH FOR JUST LOCATION=SUB

dds_sub <- DESeqDataSetFromMatrix(countData = countDatasub, colData = colDatasub, design = ~ health)

dim(dds_sub)
#[1] 13053    29

dds_sub <- dds_sub[ rowSums(counts(dds_sub)) > 100, ]
dim(dds_sub)
#[1] 12397    48  # at > 100

colData(dds_sub)$health <- factor(colData(dds_sub)$health, levels=c("H","S")) #sets that "healthy is the reference

dds_sub <- DESeq(dds_sub) 

res_sub <- results(dds_sub)
res_sub <- res_sub[order(res_sub$padj),]
head(res_sub)

summary(res_sub)
summary(res_int)
summary(res)

#################### MODEL NUMBER 2 - INTERACTIONS BETWEEN LOCATION AND HEALTH STATUS

dds_HXL <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health + location:health)

dim(dds_HXL)

dds_HXL <- dds_HXL[ rowSums(counts(dds_HXL)) > 100, ]

colData(dds_HXL)$health <- factor(colData(dds_HXL)$health, levels=c("H","S"))  #sets that "healthy is the reference

#dds_HXL <- DESeq(dds_HXL, parallel=T)
dds_HXL <- DESeq(dds_HXL, parallel=T, test="LRT", reduced = ~ health)

resultsNames(dds_HXL)
# [1]  "Intercept"           "location_sub_vs_int" "health_S_vs_H"       "locationsub.healthS"

res_HXL <- results(dds_HXL, name="locationsub.healthS")
res_HXL <- res_HXL[order(res_HXL$padj),]
head(res_HXL)

summary(res_HXL)



################################################################
# Data quality assessment by sample clustering and visualization 

par(mfrow=c(2,2))
plotMA(res, main="All samples, controlling for location", ylim=c(-2,2))
plotMA(res_int, main="Just intertidal", ylim=c(-2,2))


## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN44744_c1_g2_TRINITY_DN44744_c1_g2_i2_g.17686_m.17686", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p


## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(dds, gene="TRINITY_DN46245_c3_g3_TRINITY_DN46245_c3_g3_i2_g.21719_m.21719", intgroup=(c("health","score","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
p

############## PCA plots
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("score"))
plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))
