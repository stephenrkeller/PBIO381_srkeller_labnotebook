library(RPostgreSQL)
library(lme4)
library(pdist)
library(maps)

# For GBS-SNPs: 
GF_predict <- read.csv("~/Dropbox/GBS_Paper/data/gfGDM/offsetgfStd.csv")
GF_predictBayenv <- read.csv("~/Dropbox/GBS_Paper/data/gfGDM/offsetgfRaw.csv")

load("~/Dropbox/GBS_Paper/data/gfGDM/offsetsRandomSNPs.RData")

###########################
### Connect to PoplarDB ###
###########################

drv <- dbDriver("PostgreSQL")

con <- dbConnect(drv, host="uvm-poplar.c0l9mrq35fmv.us-east-1.rds.amazonaws.com", port=5432, user="poplar_user", password="balsamifera", dbname = "poplar")   ## creates the connection object
dbListTables(con)

#########################################################
### GF model prediction of VT garden based on FT SNPs ###
#########################################################

SPAD <- dbGetQuery(con, "SELECT pop_code,ind_code,spad,row_num,tree_num from garden_spad JOIN replicates USING(rep_code,garden_id) JOIN individual USING(ind_code) JOIN population USING(pop_code)") 

GF_SPAD <- merge(GF_predict, SPAD, by.x="pop_code", by.y="pop_code")

m1 <- lm(spad ~ VT, data=GF_SPAD)

GF_SPAD$VT2 <- GF_SPAD$VT^2
m2 <- lm(spad ~ VT + VT2, data=GF_SPAD)

predvals <- seq(0,0.12,by=0.005) 
predcounts <- predict(m1,list(VT=predvals))

pdf("~/Dropbox/GF2015/FloweringTimeSNPsManuscript/figs/GF_predict_SPAD2015.pdf", 
    width=6, height=6)
with(GF_SPAD, plot(VT, spad, ylab="Chlorophyll Fluorescence", xlab="GF Transfer Distance", main="Vermont Common Garden"))
lines(predvals, predcounts, lwd=4, col="blue")
dev.off()

######### Leaf Mass Area (LMA) #########

LMA <- dbGetQuery(con, "SELECT pop_code,ind_code,num_disks,total_weight,row_num,tree_num from _sla JOIN replicates USING(rep_code,gard_id) JOIN individual USING(ind_code) JOIN population USING(pop_code)") 
LMA$LMA <- LMA$total_weight/(31.67*LMA$num_disks) # Units mg/mm^-2

GF_LMA <- merge(GF_predict, LMA, by.x="code", by.y="pop_code")

LMA_lin <- lm(LMA ~ tDistVT, data=GF_LMA)

GF_LMA$tDistVT2 <- GF_LMA$tDistVT^2
LMA_nonlin <- lm(LMA ~ tDistVT + tDistVT2, data=GF_LMA)

predvals <- seq(0,0.19,by=0.005) 

predcounts <- predict(LMA_nonlin,list(tDistVT=predvals, tDistVT2=predvals^2))

pdf("~/Dropbox/GF2015/FloweringTimeSNPsManuscript/figs/GF_predict_LMA2015.pdf", 
    width=6, height=6)
with(GF_LMA, plot(tDistVT, LMA, ylim=c(0.02, 0.13), ylab="Leaf Mass Area (mg/mm-2)", xlab="GF Transfer Distance", main="Vermont Common Garden"))
lines(predvals, predcounts, lwd=4, col="blue")
dev.off()

######### Height Increment (HI) #########

HI <- dbGetQuery(con, "SELECT pop_code,ind_code,height_increment,row_num,
                 tree_num from garden_height JOIN replicates USING(rep_code,garden_id) JOIN individual USING(ind_code) JOIN population USING(pop_code)") 
HI$HI <- HI$height_increment

GF_HI <- merge(GF_predict, HI, by.x="pop_code", by.y="pop_code")
GF_HI_Bayenv <- merge(GF_predictBayenv, HI, by.x="pop_code", by.y="pop_code")

HI_lin_Bayenv <- lm(HI ~ VT, data=GF_HI_Bayenv)

GF_HI$VT2 <- GF_HI$VT^2
HI_nonlin <- lm(HI ~ VT + VT2, data=GF_HI)

anova(HI_lin, HI_nonlin) # n.s.

predvals <- seq(0,0.19,by=0.005) 

predcounts <- predict(HI_lin,list(VT=predvals))

pdf("~/Dropbox/GF2015/FloweringTimeSNPsManuscript/figs/GF_predict_HeightInc2015.pdf", 
    width=6, height=6)
with(GF_HI, plot(VT, HI, ylab="Height Increment (cm)", xlab="GF Transfer Distance", main="Vermont Common Garden"))
lines(predvals, predcounts, lwd=4, col="blue")
#text(.11, 110, "R-sq = 0.12",cex = 1.2)
dev.off()

# Also try first MLM to factor out garden effects and then model tDistVT as function of ind_code BLUPs
HI_MLM <- lmer(HI ~ row_num*tree_num + (1|ind_code), HI)

summary(HI_MLM)

dotplot(ranef(HI_MLM, condVar=TRUE))

ht_blups = ranef(HI_MLM) # get BLUPs from model
ht_blupsint = coef(HI_MLM) # get BLUPs from model
ht_yhat = fitted(HI_MLM) # get fitted values for levels of fixed effects
ht_resid = resid(HI_MLM) # get model residuals

ind_code <- row.names(ht_blups$ind_code)
HIblups <- as.numeric(ht_blups$ind_code[,1]) # Gets vector of BLUPs for just the individual genotypes
HIblups2 <- as.numeric(ht_blupsint$ind_code[,1]) # Gets vector of BLUPs for just the individual genotypes

htblups2 <- data.frame(ind_code, HIblups2)
htblups2$pop_code <- sub("_.+", "", htblups2$ind_code)

GF_htblups2 <- merge(GF_predict, htblups2, by.x="pop_code", by.y="pop_code")

HIblups_lin <- lm(HIblups2 ~ VT, data=GF_htblups2)

GF_htblups2$VT2 <- GF_htblups2$VT^2
HIblups_nonlin <- lm(HIblups2 ~ VT + VT2, data=GF_htblups2)

anova(HIblups_lin,HIblups_nonlin)

predvals <- seq(0,0.19,by=0.005) 

predcounts <- predict(HIblups_lin,list(VT=predvals))

pdf("~/Dropbox/GF2015/FloweringTimeSNPsManuscript/figs/GF_predict_HeightInc2015_blups_bothfits.pdf", 
    width=6, height=6)
with(GF_htblups2, plot(VT, HIblups2, ylab="Height Increment (cm)", xlab="GF Transfer Distance", main="Vermont Common Garden (Clonal BLUPs)"))
lines(predvals, predcounts, lwd=4, col="blue")
dev.off()

pdf("~/Dropbox/GF2015/FloweringTimeSNPsManuscript/figs/GF_predict_HeightInc2015_all.v.blups.pdf", 
    width=10, height=6)
par(mfrow=c(1,2))
predcounts <- predict(HI_nonlin,list(tDistVT=predvals, tDistVT2=predvals^2))
with(GF_HI, plot(tDistVT, HI, ylab="Height Increment (cm)", xlab="GF Transfer Distance", main="Vermont Common Garden (All Reps)"))
lines(predvals, predcounts, lwd=4, col="blue")

predcounts <- predict(HIblups_nonlin,list(tDistVT=predvals, tDistVT2=predvals^2))
with(GF_htblups2, plot(tDistVT, HIblups2, ylab="Height Increment (cm)", xlab="GF Transfer Distance", main="Vermont Common Garden (Clonal BLUPs)"))
lines(predvals, predcounts, lwd=4, col="blue")
dev.off()

######### Bud Flush 2016 (BF) VERMONT #########

BF <- read.table("/Volumes/kellrlab/STEVE/VT_gardenMap_budflush2016.txt", sep="\t", header=T)

GF_BF <- merge(GF_predict, BF, by.x="code", by.y="PopCode")

BF_lin <- lm(BF ~ tDistVT, data=GF_BF, na.action=na.omit)

GF_BF$tDistVT2 <- GF_BF$tDistVT^2
BF_nonlin <- lm(BF ~ tDistVT + tDistVT2, data=GF_BF, na.action=na.omit)

predvals <- seq(0,0.19,by=0.005) 

predcounts <- predict(BF_nonlin,list(tDistVT=predvals, tDistVT2=predvals^2))

with(GF_BF, plot(tDistVT, BF))
lines(predvals, predcounts, lwd=4, col="blue")

######### Phenology 2015 INDIAN HEAD #########

BF_IH <- read.csv("/Volumes/kellrlab/STEVE/scripts/QG_hsq/phenology_IH_2015.csv")
GF_BF_IH <- merge(GF_predict, BF_IH, by.x="pop_code", by.y="PopCode")

BF_IH_lin <- lm(LeafDrop ~ IH, data=GF_BF_IH, na.action=na.omit)

GF_BF_IH$tDistIH2 <- GF_BF_IH$tDistIH^2
BF_IH_nonlin <- lm(LeafDrop ~ tDistIH + tDistIH2, data=GF_BF_IH, na.action=na.omit)

predvals <- seq(0,0.21,by=0.005) 

predcounts <- predict(BF_IH_nonlin,list(tDistIH=predvals, tDistIH2=predvals^2))

with(GF_BF_IH, plot(IH, LeafDrop))
lines(predvals, predcounts, lwd=4, col="blue")

######### Bud Flush Olson (BF_Ols) INDIAN HEAD #########

BF_Ols <- read.csv("~/Dropbox/GF2015/FloweringTimeSNPsManuscript/results/gf/Phenotypes_OlsonEtAl2012.txt", sep="\t", header=T)
GF_BF_Ols <- merge(GF_predict, BF_Ols, by.x="code", by.y="pop_code")

BF_Ols_lin <- lm(IH_Budset-BudFlushIH ~ tDistIH, data=GF_BF_Ols, na.action=na.omit)

GF_BF_Ols$tDistIH2 <- GF_BF_Ols$tDistIH^2
BF_Ols_nonlin <- lm(IH_Budset-BudFlushIH ~ tDistIH + tDistIH2, data=GF_BF_Ols, na.action=na.omit)

predvals <- seq(0,0.19,by=0.005) 

predcounts <- predict(BF_Ols_nonlin,list(tDistIH=predvals, tDistIH2=predvals^2))

pdf("~/Dropbox/GF2015/FloweringTimeSNPsManuscript/figs/GF_predict_HeightGrowthDuration_IH_Olson.pdf", 
    width=6, height=6)
with(GF_BF_Ols, plot(tDistIH, IH_Budset-BudFlushIH, xlab="GF Transfer Distance", ylab="Height Growth Duration (days)", main="Indian Head Common Garden (Clonal BLUPs)"))
lines(predvals, predcounts, lwd=4, col="blue")
dev.off()

with(GF_BF_Ols, plot(tDistIH, IH_Budset-BudFlushIH, xlab="GF Transfer Distance", ylab="Height Growth Duration (days)", main="IH Common Garden"))
lines(predvals, predcounts, lwd=4, col="blue")


HI_IH <- read.csv("/Volumes/kellrlab/datashare/BalsamPoplar/UMCES/poplar__share/CommonGarden/IH/Phenology_2015/IH_phenology_2015.csv", header=T)

GF_HI_IH <- merge(GF_predict, HI_IH, by.x="pop_code", by.y="pop_code")


HI_IH_lin <- lm(height_increment ~ IH, data=GF_HI_IH, na.action=na.omit)

with(GF_HI_IH, plot(IH, height_increment, xlab="GF Transfer Distance", ylab="Height Incremenet (cm)", main="Indian Head Common Garden (Clonal BLUPs)"))

with(GF_HI_IH, plot(IH, bud_set-bud_flush, xlab="GF Transfer Distance", ylab="Height Growth Duration (days)", main="Indian Head Common Garden (Clonal BLUPs)"))

HI_IH_MLM <- lmer(height_increment ~ (1|ind_code), HI_IH)

summary(HI_MLM)

dotplot(ranef(HI_MLM, condVar=TRUE))

IHht_blups = ranef(HI_IH_MLM) # get BLUPs from model
IHht_blupsint = coef(HI_IH_MLM) # get BLUPs from model
IHht_yhat = fitted(HI_IH_MLM) # get fitted values for levels of fixed effects
IHht_resid = resid(HI_IH_MLM) # get model residuals

ind_code <- row.names(IHht_blups$ind_code)
HI_IHblups <- as.numeric(IHht_blups$ind_code[,1]) # Gets vector of BLUPs for just the individual genotypes
HIblups2 <- as.numeric(ht_blupsint$ind_code[,1]) # Gets vector of BLUPs for just the individual genotypes

HIhtblups2 <- data.frame(ind_code, HI_IHblups)
HIhtblups2$pop_code <- sub("_.+", "", HIhtblups2$ind_code)

IH_GF_htblups2 <- merge(GF_predict, HIhtblups2, by.x="pop_code", by.y="pop_code")

HIblups_lin <- lm(HI_IHblups ~ IH, data=IH_GF_htblups2)

IH_GF_htblups2$IH2 <- IH_GF_htblups2$IH^2
HIblups_nonlin <- lm(HI_IHblups ~  IH + IH2, data=IH_GF_htblups2)

# Try looping over all the 999 random SNP draws from MCF and do regression:

# Try one: nope
#models <- lapply(offsetgfRand, summary(lm, formula= VT~IH))

# Try two:
est=numeric(length=999)
for (i in 1:length(offsetgfRand)) {
  df <- merge(HI, offsetgfRand[[i]], by.x="pop_code", by.y="pop_code")
  mod <- lm(df$height_increment ~ df$VT)
  est[i] <- mod$coefficients[2]
}

hist(unlist(est), breaks=50, xlim=c(-1200,100), col="gray")
abline(v=HI_lin$coefficients[2], col="red", lwd=2)





