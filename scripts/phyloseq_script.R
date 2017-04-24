library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
library("ggplot2")
theme_set(theme_bw())

#Import the OTU table
otutable <- import_biom(BIOMfilename = 'otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom', 
                        treefilename = 'rep_set_no_chimeras.tre', 
                        parseFunction = parse_taxonomy_greengenes)

#The warnings are ok. There is 1 warning for every OTU that doesn't have a taxonomy assignment

#Import the mapping file
mapping <- import_qiime_sample_data(mapfilename = 'R_map.txt')

#Merge the mapping file to the OTU table
phylo <- merge_phyloseq(otutable, mapping)

#Check to make sure the imports worked
phylo

###############################################################################
####Testing for differentially expressed OTUs between Sick and Healthy samples
####while controlling for the repeated measures on each individual
###############################################################################

##For now we're only going to look at differences between Sick and Healthy samples and 
##remove the Dead samples from the analysis
phylo_subset = subset_samples(phylo, Phenotype != "Dead")

##We want to numbers of each individual to be a factor. Right now it's an integer so we have to change that.
class(sample_data(phylo_subset)$individual)
sample_data(phylo_subset)$individual<-factor(sample_data(phylo_subset)$individual)

##Phyloseq's wrapper to get OTU data into DESeq
pheno <- phyloseq_to_deseq2(phylo_subset, ~ individual + Phenotype)

##Run DESeq. This command takes a bit of time
pheno_deseq_test <- DESeq(pheno, test="Wald")

##Get results from DESeq
pheno_results <- results(pheno_deseq_test)

##Have a look at the summary
summary(pheno_results)
head(pheno_results)

##Make table with OTUs padj<0.05
alpha <- 0.05
pheno_sigtab <- pheno_results[which(pheno_results$padj < alpha), ]

##Add taxa info to that table
pheno_sigtab <- cbind(as(pheno_sigtab, "data.frame"), as(tax_table(phylo)[rownames(pheno_sigtab), ], "matrix"))
tail(pheno_sigtab)

##Save the table to your desktop
write.table(pheno_sigtab, "DE_OTU_sickk_vs_healthy.txt", sep="\t")

##plot taxa of all otus padj<0.05 using ggplots2
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}

x <- tapply(pheno_sigtab$log2FoldChange, pheno_sigtab$Phylum, function(x) max(x))
x <- sort(x, TRUE)
pheno_sigtab$Phylum = factor(as.character(pheno_sigtab$Phylum), levels=names(x))
# Family order
x <- tapply(pheno_sigtab$log2FoldChange, pheno_sigtab$Family, function(x) max(x))
x <- sort(x, TRUE)
pheno_sigtab$Family = factor(as.character(pheno_sigtab$Family), levels=names(x))
ggplot(pheno_sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

##plot counts across days for individual otus
title <- "New.ReferenceOTU2814"
data <- plotCounts(pheno_deseq_test, "New.ReferenceOTU2814" , 
                   intgroup=c("Day","Phenotype"), returnData=TRUE)
ggplot(data, aes(x=Day, y=count, color=Phenotype, group=Phenotype)) + 
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() + ggtitle(title)


######################################################################################
###Testing for differentially expressed between samples with different symptom numbers
######################################################################################
sample_data(phylo)$individual <- factor(sample_data(phylo)$individual)
sample_data(phylo)$Pheno_num <- factor(sample_data(phylo)$Pheno_num)

##Phyloseq's wrapper to get OTU data into DESeq
pheno_num <- phyloseq_to_deseq2(phylo, ~ individual + Pheno_num)

##Run DESeq. This command takes a bit of time
pheno_num_deseq_test <- DESeq(pheno_num, test="Wald")

##Get results from DESeq
pheno_num_results <- results(pheno_num_deseq_test)
head(pheno_num_results)

##Compare Healthy to S_1
pheno_num_res_0_1<- results(pheno_num_deseq_test, contrast=c("Pheno_num","0","1"))
head(pheno_num_res_0_1)
summary(pheno_num_res_0_1)

##Compare Healthy to S_2
pheno_num_res_0_2 <- results(pheno_num_deseq_test, contrast=c("Pheno_num","0","2"))
head(pheno_num_res_0_2)
summary(pheno_num_res_0_2)
#Make a table of significant resutls and add the taxonomic information
alpha <- 0.05
results_table_0_2 <- pheno_num_res_0_2[which(pheno_num_res_0_2$padj < alpha), ]
results_table_0_2 <- cbind(as(results_table_0_2, "data.frame"), as(tax_table(phylo)[rownames(results_table_0_2), ], "matrix"))
head(results_table_0_2)
#Save it to your desktop if you want

##Compare Healthy to S_3
pheno_num_res_0_3 <- results(pheno_num_deseq_test, contrast=c("Pheno_num","0","3"))
head(pheno_num_res_0_3)
summary(pheno_num_res)
#Make the table of significant results as we did above

##Compare Healthy to S_4
pheno_num_res_0_4<- results(pheno_num_deseq_test, contrast=c("Pheno_num","0","4"))
head(pheno_num_res_0_4)
summary(pheno_num_res_0_4)

##Compare Healthy to S_5
pheno_num_res_0_5<- results(pheno_num_deseq_test, contrast=c("Pheno_num","0","5"))
head(pheno_num_res_0_5)
summary(pheno_num_res)

##Plot abundances of taxa in different pheno numbers
##First we must rarefy
set.seed(28132)
phyloR = rarefy_even_depth(phylo, sample.size = 18000)

#Then we check to make sure we rarefied corectly.
title = "Sum of reads for each sample, phyloR"
plot(sort(sample_sums(phyloR), TRUE), type = "h", main = title, ylab = "reads", 
    ylim = c(0, 20000))


##We merge on basis of pheno number so we can make relative abundance graph
phyloRm = merge_samples(phyloR, "Pheno_num")

phyloRt = transform_sample_counts(phyloRm, function(x) 100 * x/sum(x))

p = plot_bar(phyloRt, "Pheno_num", fill="Order")
p + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

##Merge on basis of phenotype
phyloRmD = merge_samples(phyloR, "Final_phenotype")
phyloRtD = transform_sample_counts(phyloRmD, function(x) 100 * x/sum(x))

##Make the plot
p = plot_bar(phyloRtD , "Final_phenotype", fill="Family")
p + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")



