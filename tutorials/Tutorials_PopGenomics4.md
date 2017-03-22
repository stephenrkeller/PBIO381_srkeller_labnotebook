P/BIO381 Tutorials

## Population Genomics 4: Population Structure with PCA and ADMIXTURE

### March 22, 2017

Our next goal is to look for the presence of population structure in our sample of sea stars. Recall that these animals were all collected from the same general geographic area, and the dispersal ability of sea star gametes and juvelines is pretty impressive. So, we don't necessarily expect to find a lot of structure, but one nevers knows without checking...

We'll take 2 different approaches to test if there is any population structure present in our sample: 

1. A Principal Component Analysis (PCA) on the SNPs to see if they group by sampling locality

2. We'll use the ADMIXTURE program to cluster genotypes into  *K*  groups, in which we'll vary *K* from 1 - 10. 

   ​

   Keep in mind, both analyses are naive with regard to the actual sampling locality of individuals, so they provide a relatively unbiased way of determining if there are actually >1 genetically distinct groups represented in the data.

------------------------------

## PCA on SNP genotypes:##

Principal Components Analysis (PCA) is a powerful multivariate technique to reduce the dimensionality of large SNP datasets into a few synthetic axes (PC's) that describe the major structure present in the data. We'll do this in **R** using the *adegent* package ([adegenet manual available here](https://cran.r-project.org/web/packages/adegenet/adegenet.pdf)).



* Transfer your filtered vcf.gz file (SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz) from the server to your local machine. You know the drill…use Fetch, WinScp, or scp at the command-line.
* Also transfer the metadata on ssw locality and disease, found here:

```
/data/project_data/snps/reads2snps/ssw_healthloc.txt
```



* Open **R**, paste the following into an R script, and work through it:

```R
# Set your working directory to where you downloaded your results files:
setwd("~/github/PBIO381_srkeller_labnotebook/data/SNP_data/")

list.files() # Do you see your downloaded files there? If not, double check to make sure you've set your working directory to the right spot

# We'll need to install 2 packages to work with the SNP data:
install.packages("vcfR") # reads in vcf files and proides tools for file conversion 
install.packages("adegenet") # pop-genetics package with some handy routines, including PCA and other multivariate methods (DAPC)

# ...and load the libraries
library(adegenet)
library(vcfR)

#Read the vcf SNP data into R
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")

# The adegenet package uses a highly efficient way of storing large SNP datasets in R called a "genlight" object. The following function creates a genlight object from your vcf:
gl1 <- vcfR2genlight(vcf1)
print(gl1) # Looks good! Right # of SNPs and individuals!

# For info, try:
gl1$ind.names
gl1$loc.names[1:10]

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort it by Individual ID

# Confirm the ID's are ordered the same in gl1 and ssw_meta:
gl1$ind.names
ssw_meta$Individual

gl1$pop <- ssw_meta$Locality # assign locality info
gl1$other <- as.list(ssw_meta$Trajectory) # assign disease status

# Now, let's compute the PCA on the SNP genotypes and plot it:
PCA1 <- glPca(gl1, nf=4) 
PCA1 # prints summary
scatter(PCA1, label=gl1$pop) # plots PCA scores and labels by locale
loadingplot(pca1) # which SNPs load most strongly on the 1st PC axis?
```

------------------------------------------

## ADMIXTURE analysis ##

A second way of estimating population structure besides PCA is to use genotypic clustering algorithms. These include the familiar program STRUCTURE, as well as many others that have sprung up like it. All share the common feature of using multi-locus genetic data to estimate:

* (i) the number of clusters present, and 
* (ii) each individual's proportion of genetic ancestry in these clusters

With large population genomic datasets, STRUCTURE would take a prohibitively long time to run. Thus, analyzing thousands to millions of SNPs requires computationally efficient approaches to the clustering problem. A good option is the maximum-likelihood program ADMIXTURE by John Novembre's lab.

For reference, here is the source page for information on [ADMIXTURE](https://www.genetics.ucla.edu/software/admixture/).

And as with any good software, there is also a well annotated [manual](https://www.genetics.ucla.edu/software/admixture/admixture-manual.pdf) available.

ADMIXTURE introduces a user-defined number of groups or clusters (known as K) and uses maximum likelihood to estimate allele frequencies in each cluster, and assign each individual ancestry (Q) to one or more of these clusters. 

From a practical standpoint, ADMIXTURE is pretty easy to run. At a minimum, it just needs an input file and the requested level of K to investigate. Unfortunately, getting the data formatted properly for input is a real pain for non-model organisms. 

I used the program [pgdspider](http://www.cmpg.unibe.ch/software/PGDSpider/) to convert from our vcf files to .geno. The ready to go file is located on our server here:

```
/data/project_data/snps/reads2snps/SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno
```

In the same path, you should also see a bash script file called:  `ADMIX.sh`

Use **cp** to copy the .geno and ADMIX.sh files to your *home directory on the server*, then **cd** there and confirm the files are present.

From within your home directory, open the ADMIX.sh script in vim. Let's walk through what each step is doing:

```bash
#!/bin/bash

# Run ADMIXTURE to determine the number of genetic clusters in the SNP data, 
# and the ancestry proportions of each individual

# Here's the utility of 'for loops'...

for K in {1..10}

do

admixture -C 0.000001 --cv ./SSW_all_biallelic.MAF0.02.Miss1.0.recode.vcf.geno $K \
| tee log${K}.out

done

grep "CV" log*.out >chooseK.txt
```



When you're ready to go, exit vim to return to the command line, and execute the script.

```bash
$ bash ADMIX.sh
```



The cross-validation procedure in ADMIXTURE breaks the samples into 5 equally sized chunks. It then masks each chunk in turn, trains the model to estimate the allele frequencies and ancestry assignments on the unmasked data, and then attempts to predict the genotype values for the masked individuals. 

**If the model is good (and there's true structure in the data), then the best value of K is the one that will *minimize* the cross-validation (CV) error.**

![ADMIXTURE CV](https://www.researchgate.net/profile/Jason_Hodgson/publication/263579532/figure/download/fig3/AS:392426666643462@1470573216485/Figure-S1-Plot-of-ADMIXTURE-cross-validation-error-from-K2-through-K6-We-chose-K3-to.png)The CV values for our runs are stored in the output file "chooseK.txt"

Print the contents of this file to your screen:

```bash
$ cat chooseK.txt
```

* What level of K is the CV the lowest? 
* What does this say about the presence of genetic structure in our SSW data?



