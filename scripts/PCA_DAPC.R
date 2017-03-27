setwd("~/github/PBIO381_srkeller_labnotebook/data/SNP_data")

list.files()

#install.packages("vcfR")
#install.packages("adegenet")

library(vcfR)
library(adegenet)

vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
ssw_meta <- read.table("ssw_healthloc.txt", header=T)
ssw_meta <- ssw_meta[order(ssw_meta$Individual),]

gl1 <- vcfR2genlight(vcf1)
gl1$pop <- ssw_meta$Location
gl1$other <- as.list(ssw_meta)

glPlot(gl1, posi="bottomleft") #plots heatmap of samples x loci

pca1 <- glPca(gl1, nf=4) # Perform PCA on SNP genotypes
pca1 # gives summary

#scatter(pca1, posi="bottomleft", label=gl1$pop) # funky plot with block labels
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=gl1$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(gl1$pop), 
       pch=20, 
       col=c("black", "red"))

plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=as.factor(gl1$other$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (5317 SNPs)")
legend("topleft", 
       legend=unique(gl1$other$Trajectory), 
       pch=20, 
       col=as.factor(unique(gl1$other$Trajectory)))


loadingplot(abs(pca1$loadings[,1]))  
loadingplot(abs(pca1$loadings[,1]), threshold=quantile(abs(pca1$loadings), 0.999))
gl1$loc.names[which(abs(pca1$loadings)>0.1)]

vcf2 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss1.0.recode.vcf")

gl2 <- vcfR2genlight(vcf2)
gl2$pop <- ssw_meta$Location
gl2$other <- as.list(ssw_meta)

glPlot(gl2, posi="bottomleft")

pca2 <- glPca(gl2, nf=4)
pca2 # prints summary

plot(pca2$scores[,1], pca2$scores[,2], 
     cex=2, pch=20, col=as.factor(gl1$pop), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=0%; 1494 SNPs)")
legend("topleft", 
       legend=unique(gl2$pop), 
       pch=20, 
       col=as.factor(unique(gl1$pop)))

plot(pca2$scores[,1], pca2$scores[,2], 
     cex=2, pch=20, col=as.factor(gl2$other$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (1494 SNPs)")
legend("topleft", 
       legend=unique(gl2$other$Trajectory), 
       pch=20, 
       col=as.factor(unique(gl2$other$Trajectory)))

gl2$other$dispop <- paste(gl2$pop,gl2$other$Trajectory)

plot(pca2$scores[,1], pca2$scores[,2], 
     cex=2, pch=(15+as.numeric(gl2$pop)), 
     col=as.factor(gl1$other$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (1494 SNPs)")

legend("topleft", 
       legend=levels(as.factor(gl2$other$dispop)), 
       pch=c(16, 16, 16, 17, 17, 17), 
       col=as.factor(gl2$other$Trajectory))

loadingplot(abs(pca2$loadings[,1]), threshold=quantile(abs(pca2$loadings), 0.999))
loadingplot(abs(pca2$loadings[,1]), fac=gl2$chromosome, byfac=T, threshold=0.5)



# DAPC
# On tidal groups
d1 <- dapc(gl2, pop=gl1$pop, n.pca=8, n.da=4,
           var.loadings=T, pca.info=T,
           pca.select=c("nbEig", "percVar"))

scatter.dapc(d1, grp=gl1$pop, legend=T)
compoplot(d1)
assignplot(d1)

# On disease groups
d2 <- dapc(gl2, pop=gl1$other$Trajectory, n.pca=8, n.da=4,
     var.loadings=T, pca.info=T,
     pca.select=c("nbEig", "percVar"))

scatter.dapc(d2, grp=gl1$other$Trajectory, legend=T)
compoplot(d2)
#assignplot(d2)

loadingplot(abs(d2$var.load, lab.jitter=1, threshold=quantile(abs(d2$var.load), probs=0.999))



