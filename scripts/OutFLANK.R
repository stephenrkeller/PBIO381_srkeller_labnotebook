install.packages("devtools")
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
install_github("whitlock/OutFLANK")

library(OutFLANK)
library(vcfR)

setwd("~/github/PBIO381_srkeller_labnotebook/data/SNP_data")

ssw.geno <- read.fwf("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno", width=rep(1,24))
ssw_meta <- read.table("ssw_healthloc.txt", header=T)
ssw_meta <- ssw_meta[order(ssw_meta$Individual),]

ssw_meta$Disease <- factor(rep(NA, length(ssw_meta$Trajectory)), levels=c("Healthy", "Sick") )   
ssw_meta$Disease[ssw_meta$Trajectory %in% c("HS", "SS")] <- "Sick"
ssw_meta$Disease[ssw_meta$Trajectory %in% "HH"] <- "Healthy"

ssw_meta$Trajectory[which(ssw_meta$Trajectory=="MM")] = NA

ssw.geno.trans <-t(ssw.geno)

OF_SNPs <- MakeDiploidFSTMat(ssw.geno.trans, locusNames=seq(1, 5317, by=1), popNames=ssw_meta$Disease)

OF_out <- OutFLANK(FstDataFrame=OF_SNPs, LeftTrimFraction=0.1, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=3, qthreshold=0.1)

OutFLANKResultsPlotter(OF_out, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, Zoom=F, RightZoomFraction=0.05, titletext=NULL)

outliers <- which(OF_out$results$OutlierFlag=="TRUE")
outliers

vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
head(getFIX(vcf1))
vcfann <- as.data.frame(getFIX(vcf1))
vcfann[outliers,]
