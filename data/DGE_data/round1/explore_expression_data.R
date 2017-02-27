

setwd("~/Dropbox/1_Research/SSW/DGE")

countsTable <- read.delim('countstatsummary.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
attach(countsTable)
head(countsTable)


plot(NumMultiAligned,NumSingleAligned, xlim=c(0,4000000), ylim=c(0,4000000))
abline(1,1)

countsTable <- read.delim('countssummarystats_beetle.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
attach(countsTable)
head(countsTable)
max(NumContigsMatched)
