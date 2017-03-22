setwd("~/github/PBIO381_srkeller_labnotebook/data/SNP_data")

list.files()

install.packages("vcfR")
install.packages("adegenet")

library(adegenet)
library(vcfR)
library(SNPRelate)

vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
ssw_meta <- read.table("../ssw_healthloc.txt", header=T)
ssw_meta <- ssw_meta[order(ssw_meta$Individual),]

gl1 <- vcfR2genlight(vcf1)
gl1$pop <- ssw_meta$Location
gl1$other <- as.list(ssw_meta)

#glPlot(gl1, posi="topleft")

pca1 <- glPca(gl1, nf=4)
pca1 # prints summary

scatter(pca1, posi="bottomleft", label=gl1$pop)
loadingplot(pca1)

# DAPC
d1 <- dapc(gl1, pop=gl1$pop, n.pca=6, n.da=4,
     scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
     pca.select=c("nbEig", "percVar"))
scatter.dapc(d1, grp=gl1$other$pop, scree.da=F, scree.pca=F)

loadingplot(d1$var.contr, lab.jitter=1)
head(d1$var.contr)

# Compare to PCA in SNPRelate package





