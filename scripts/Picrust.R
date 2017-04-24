library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
library("ggplot2")
theme_set(theme_bw())
library('biom')

x = read_biom("metagenome_predictions.L3.biom")
otumat = as(biom_data(x), "matrix")
OTU = otu_table(otumat, taxa_are_rows=TRUE)


mapping <- import_qiime_sample_data(mapfilename = 'R_map.txt')

phylo <- merge_phyloseq(OTU, mapping)
phylo

###############################################################################
###Test for DE KO terms between individuals that got sick and those that didn't
###############################################################################

final_pheno = phyloseq_to_deseq2(phylo, ~ Final_phenotype)
final_pheno_results = DESeq(final_pheno, test="Wald")
final_pheno_res = results(final_pheno_results)
summary(final_pheno_res)
head(final_pheno_res)

alpha = 0.05
final_pheno_sigtab = final_pheno_res[which(final_pheno_res$padj < alpha), ]
final_pheno_sigtab= cbind(as(final_pheno_sigtab, "data.frame"), as(tax_table(phylo)[rownames(final_pheno_sigtab), ], "matrix"))
head(final_pheno_sigtab)
final_pheno_sigtab
write.table(final_pheno_sigtab, "Final_pheno_L3.txt", sep="\t")