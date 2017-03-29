setwd("~/github/PBIO381_srkeller_labnotebook/results")

# Import the ADMIXTURE Q matrices
K1Q <- read.table("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.1.Q")
K2Q <- read.table("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.2.Q")
K3Q <- read.table("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.3.Q")

# Get the SSW meta-data
ssw_meta <- read.table("~/github/PBIO381_srkeller_labnotebook/data/SNP_data/ssw_healthloc.txt", header=)
  
# Set up the plotting coditioms for a multi-panel plot (4 rows, 1 column)
par(mfrow=c(3,1))

# Make the barplots for K=1-3
barplot(t(as.matrix(K1Q)), 
        col=rainbow(2),
        names.arg=ssw_meta$Location, 
        cex.names=0.75, 
        xlab="Individual", ylab="Ancestry", 
        border=NA)
barplot(t(as.matrix(K2Q)), 
        col=rainbow(2),
        names.arg=ssw_meta$Location, 
        cex.names=0.75, 
        xlab="Individual", ylab="Ancestry", 
        border=NA)
barplot(t(as.matrix(K3Q)), 
        col=rainbow(3),
        names.arg=ssw_meta$Location, 
        cex.names=0.75, 
        xlab="Individual", ylab="Ancestry", 
        border=NA)
