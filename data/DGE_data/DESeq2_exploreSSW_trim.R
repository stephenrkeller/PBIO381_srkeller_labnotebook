source("http://bioconductor.org/workflows.R")
workflowInstall("rnaseqGene")

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library(DESeq2)
library(ggplot2)

setwd("~/github/PBIO381_srkeller_labnotebook/data/DGE_data")

countsTable <- read.delim('countsdata_trim.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

#################### Build dataset, model, and run analyses

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ day + location + health)
# In this typical model, the "sex effect" represents the overall effect controlling for differences due to population and devstage. page 26-27 manual DESeq2.pdf
# The last term in the model is what is tested.  In this case sex.
#  This is not the same as an interaction.


dim(dds)
[1] 26550    65

dds <- dds[ rowSums(counts(dds)) > 100, ] # This subsets the data for only genes that have a summed total of 100 counts or more
dim(dds)
[1] 13334    65  # at > 100; little more than an average of 10 reads per sample for the 93 samples

colSums(counts(dds))
hist(colSums(counts(dds)), breaks = 80, xlim=c(0,max(colSums(counts(dds)))))


colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S"))

dds <- DESeq(dds)  # this step takes a loooong time ~4 minutes with the trimmed data set
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 3308 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

save(dds, file="dds.trim.Robject")

res <- results(dds)
res <- res[order(res$padj),]
head(res)
# log2 fold change (MAP): health S vs H 
# Wald test p-value: health S vs H 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE
# <numeric>      <numeric> <numeric>
# TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i1_g.23138_m.23138 1950.0719       2.488783 0.4311875
# TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110  902.2693       2.475891 0.4599085
# TRINITY_DN43359_c0_g1_TRINITY_DN43359_c0_g1_i1_g.14658_m.14658  889.9707       1.163219 0.2482335
# TRINITY_DN47215_c1_g4_TRINITY_DN47215_c1_g4_i3_g.25054_m.25054  774.1126       1.723917 0.3650258
# TRINITY_DN47215_c0_g1_TRINITY_DN47215_c0_g1_i5_g.25051_m.25051  911.7634       1.586693 0.3431307
# TRINITY_DN45416_c4_g2_TRINITY_DN45416_c4_g2_i3_g.19333_m.19333 1629.8753       1.775765 0.3873817
# stat       pvalue         padj
# <numeric>    <numeric>    <numeric>
# TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i1_g.23138_m.23138  5.771927 7.837024e-09 8.691260e-06
# TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110  5.383443 7.307426e-08 4.051968e-05
# TRINITY_DN43359_c0_g1_TRINITY_DN43359_c0_g1_i1_g.14658_m.14658  4.685987 2.786136e-06 7.724563e-04
# TRINITY_DN47215_c1_g4_TRINITY_DN47215_c1_g4_i3_g.25054_m.25054  4.722727 2.327027e-06 7.724563e-04
# TRINITY_DN47215_c0_g1_TRINITY_DN47215_c0_g1_i5_g.25051_m.25051  4.624166 3.761091e-06 8.342101e-04
# TRINITY_DN45416_c4_g2_TRINITY_DN45416_c4_g2_i3_g.19333_m.19333  4.584018 4.561241e-06 8.430694e-04

summary(res)
# out of 13334 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 50, 0.37% 
# LFC < 0 (down)   : 8, 0.06% 
# outliers [1]     : 539, 4% 
# low counts [2]   : 11686, 88% 
# (mean count < 80)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

plotMA(res, main="DESeq2", ylim=c(-2,2))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN46709_c0_g1_TRINITY_DN46709_c0_g1_i1_g.23138_m.23138", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = date)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + scale_y_log10(breaks=c(25,100,1000)) + ylim(0,2500)
p

# 2.2 Data quality assessment by sample clustering and visualization 

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))

# rld <- rlog(dds, blind=FALSE) # this takes too long with such a large data set!
# plotPCA(rld, intgroup=c("health","date"))


#to save the plot as a pdf
pdf(file="PCA_1v2_allgenes.pdf", height=5.5, width=5.5)
plotPCA(rld, intgroup=c("health","day","location"))
dev.off()

pdf(file="PCA_1v2.pdf", height=5.5, width=5.5)
data <- plotPCA(rld, intgroup=c("health","day","location"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=sex, shape=devstage)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()




#########
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE)
pheatmap(assay(vsd)[select,], show_rownames=FALSE)


################ testing for interactions between factors

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ health + location + health:location)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)


dds <- DESeq(dds, parallel=T)


resultsNames(dds)
# [1] "Intercept"           "health_S_vs_H"       "location_sub_vs_int" "healthS.locationsub"

res <- results(dds, name="healthS.locationsub")
res <- res[order(res$padj),]
head(res)
# log2 fold change (MLE): healthS.locationsub 
# Wald test p-value: healthS.locationsub 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE
# <numeric>      <numeric> <numeric>
#   TRINITY_DN44444_c9_g1_TRINITY_DN44444_c9_g1_i4_g.17034_m.17034 126.60114     -28.469401  2.918454
# TRINITY_DN2191_c0_g1_TRINITY_DN2191_c0_g1_i1_g.232_m.232        42.12038     -22.438459  3.335695
# TRINITY_DN41408_c0_g3_TRINITY_DN41408_c0_g3_i1_g.11127_m.11127  23.90548     -21.338403  3.166769
# TRINITY_DN46124_c1_g2_TRINITY_DN46124_c1_g2_i6_g.21322_m.21322 414.33733      -8.743311  1.367619
# TRINITY_DN47096_c0_g1_TRINITY_DN47096_c0_g1_i7_g.24574_m.24574 106.46107      -9.946944  1.616331
# TRINITY_DN35881_c0_g1_TRINITY_DN35881_c0_g1_i1_g.5955_m.5955    28.47358     -21.035939  3.441188
# stat       pvalue         padj
# <numeric>    <numeric>    <numeric>
#   TRINITY_DN44444_c9_g1_TRINITY_DN44444_c9_g1_i4_g.17034_m.17034 -9.754959 1.756711e-22 1.468610e-18
# TRINITY_DN2191_c0_g1_TRINITY_DN2191_c0_g1_i1_g.232_m.232       -6.726773 1.734674e-11 4.833957e-08
# TRINITY_DN41408_c0_g3_TRINITY_DN41408_c0_g3_i1_g.11127_m.11127 -6.738225 1.603329e-11 4.833957e-08
# TRINITY_DN46124_c1_g2_TRINITY_DN46124_c1_g2_i6_g.21322_m.21322 -6.393089 1.625680e-10 3.397670e-07
# TRINITY_DN47096_c0_g1_TRINITY_DN47096_c0_g1_i7_g.24574_m.24574 -6.154026 7.554029e-10 1.263034e-06
# TRINITY_DN35881_c0_g1_TRINITY_DN35881_c0_g1_i1_g.5955_m.5955   -6.112987 9.778343e-10 1.362449e-06

summary(res)
# out of 13321 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 2, 0.015% 
# LFC < 0 (down)   : 111, 0.83% 
# outliers [1]     : 463, 3.5% 
# low counts [2]   : 4511, 34% 
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


d <- plotCounts(dds, gene="TRINITY_DN46124_c1_g2_TRINITY_DN46124_c1_g2_i6_g.21322_m.21322", intgroup=(c("location","health","day")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location, colour = day, fill=location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100,4000)) + scale_x_discrete(limits=c("H","S")) + scale_shape_manual(values=c(21,24))
p

p <- ggplot(d, aes(x= health, y=count, shape = location, fill=location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100,4000)) + scale_x_discrete(limits=c("H","S")) + scale_shape_manual(values=c(21,24))
p

ggsave("dot_plot-TRINITY_DN46124_c1_g2_TRINITY_DN46124_c1_g2_i6_g.21322_m.21322.png", p, width=8, height=4, dpi=300)

