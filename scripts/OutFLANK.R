install.packages("devtools")
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
install_github("whitlock/OutFLANK")

library(OutFLANK)

ssw.geno <- read.fwf("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno", width=rep(1,24))
ssw_meta <- read.table("ssw_healthloc.txt", header=T)

ssw.geno.trans <-t(ssw.geno)

OF_SNPs <- MakeDiploidFSTMat(ssw.geno.trans, locusNames=seq(1, 5317, by=1), popNames=ssw_meta$Location)

OF_out <- OutFLANK(FstDataFrame=OF_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.01, Hmin=0.01, NumberOfSamples=2, qthreshold=0.1)

OutFLANKResultsPlotter(OF_out, withOutliers=T, NoCorr=T, Hmin=0.01, binwidth=0.005, Zoom=F, RightZoomFraction=0.05, titletext=NULL)

outliers <- which(OF_out$results$OutlierFlag=="TRUE")
outliers


