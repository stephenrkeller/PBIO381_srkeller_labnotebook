library(pcadapt)
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)
setwd("~/github/PBIO381_srkeller_labnotebook/data/SNP_data")

ssw_snps <- read.pcadapt("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf", type="vcfR", allele.sep="|")
ssw_meta <- read.table("ssw_healthloc.txt", header=T)

ssw_pca <- pcadapt(ssw_snps, K=10, min.maf=0.05)

plot(ssw_pca, option="screeplot")
plot(ssw_pca, option="scores", pop=ssw_meta$Trajectory)
plot(ssw_pca)

ssw_K1 <- pcadapt(ssw_snps, K=1, min.maf=0.05)

summary(ssw_K1)
plot(ssw_K1, option="manhattan")
plot(ssw_K1, option="qqplot", threshold=0.1)
plot(ssw_K1, option="stat.distribution")

qval <- qvalue(ssw_K1$pvalues)$qvalues
alpha = 0.1
print(outliers <- which(qval<alpha))



