# Set your working directory to where you downloaded your results files:
setwd("~/github/PBIO381_srkeller_labnotebook/data/SNP_data/")

list.files() # Do you see your downloaded files there? If not, double check to make sure you've set your working directory to the right spot

# We'll need to install 2 packages to work with the SNP data:
install.packages("vcfR") # reads in vcf files and proides tools for file conversion 
install.packages("adegenet") # pop-genetics package with some handy routines, including PCA and other multivariate methods (DAPC)

library(vcfR)
library(adegenet)

vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")

gl1 <- vcfR2genlight(vcf1)
print(gl1)
gl1$ind.names
gl1$loc.names[1:10]
gl1$chromosome[1:3]

ssw_meta <- read.table("ssw_healthloc.txt", header=T) 
ssw_meta <- ssw_meta[order(ssw_meta$Individual),]

gl1$pop <- ssw_meta$Location
gl1$other <- as.list(ssw_meta$Trajectory) 

glPlot(gl1, posi="bottomleft")

pca1 <- glPca(gl1, nf=4, parallel = F) # If a PC user, then add parallel = F
pca1 

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
     cex=2, pch=20, col=as.factor(unlist(gl1$other)), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(as.factor(unlist(gl1$other))), 
       pch=20, 
       col=as.factor(unique(unlist(gl1$other))))

loadingplot(abs(pca1$loadings[,1]),
            threshold=quantile(abs(pca1$loadings), 0.999))

threshold=quantile(abs(pca1$loadings), 0.999)

gl1$loc.names[which(abs(pca1$loadings)>threshold)]

disease.dapc <- dapc(gl1, pop=as.factor(unlist(gl1$other)), n.pca=8, n.da=3,
                     var.loadings=T, pca.info=T)


scatter.dapc(disease.dapc, grp=as.factor(unlist(gl1$other)), legend=T)

loadingplot(abs(disease.dapc$var.load), 
            lab.jitter=1, 
            threshold=quantile(abs(disease.dapc$var.load), probs=0.999))

compoplot(disease.dapc)
?dapc
