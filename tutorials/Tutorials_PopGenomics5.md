P/BIO381 Tutorials

## Population Genomics 5: Testing for selective sweeps associated with disease status

### April 3, 2017 

As we've discussed in class, positive directional selection will increase the frequency of a beneficial SNP allele within a population. This often causes variation at nearby linked SNP loci to also increase in frequency along with the beneficial SNP due to LD — this is what we call a "selective sweep". There are several statistical "footprints" that selective sweeps are predicted to leave in SNP frequency data:

* reduced diversity (pi)
* increased linkage disequilibrium (LD) around the selected site
* high frequency of rare variants due to recent mutations following a sweep (skewed site frequency spectrum; negative Tajima's D)

So far, we've been thinking about detecting sweeps within a single population. However, sweeps that are *happening differently across popuations* are expected to lead to an additional signal of high differentiation between groups (i.e., high Fst). This is the basis of our next analysis approach aimed at detecting evidence of selection operating differently between healthy and sick sea stars.

------------------------------

## Challenges for finding Fst outliers:##

Fst measures the degree of allele frequency differences between 2 or more groups of individuals. In most cases, these groups represent different spatially distinct populations. However, population demographic history (drift, bottlenecks, range expansions) can have a profound effect on the distribution of Fst, and often leads to increases in the number of high Fst loci, leading to false positives in scans for selective sweeps. Therefore, we need to be careful to separate out loci with high Fst due to neutral demographic history from those with high Fst due to natural selection. 

On such approach that we'll use here was [recently published by Mike Whitlock and Katie Lotterhos](http://www.journals.uchicago.edu/doi/abs/10.1086/682949), called 'OutFLANK'. The overall goal is to use the expected neutral distribution of Fst between *n* groups based on established population genetic theory to identify observed loci that exceed the upper neutral distribution — these are our candidates for local selection. However, the basic population genetic model used to generate Fst  (the "Island Model") will fail under most realistic demographic scenarios that don't strictly follow the Island Model's assumptions of constant equal population sizes with equal and symmetric gene flow between groups. 

OutFLANK's approach is to calculate the theoretical distribution of Fst that provides the best fit to *just the neutral loci in our sample*, given the (unknown) demographic history of our species. How does it do that? First, it trims out loci in the upper and lower tails of the distribution (these are most likely to be under selection). Then, it uses the remaining loci to fit a neutral distribution of Fst based on popgen theory. It then checks to see if this neutral distribution provides a good fit to the trimmed distribution. If it does, then it identifies outlier loci beyond the upper tail of the neutral distribution as candidates for local adaptation. 

![](http://www.nature.com/scitable/content/ne0000/ne0000/ne0000/ne0000/15836493/f1_nosil.jpg)

The manual for the OutFLANK **R** package can be found [here](https://github.com/whitlock/OutFLANK/blob/master/OutFLANK%20readme.pdf).

For our analysis, our "populations" will consist of 3 groups of samples representing the disease status of sea stars: HH, HS, or SS (we'll exclude the ambiguous MM individuals here). 

The data format to read into OutFLANK is a bit funny. It is very similar to the .geno file we generated for our ADMIXTURE analysis, where SNP genotypes are scored as 0,1,2 (missing =9). The difference is that OutFLANK expects this data matrix inverted, with individuals in rows and SNPs in columns. For us, this will be a 22 x 5317 matrix (22 because we're not including the MM's).



OK, let's get started:

* Download your "SSW_all_biallelic.MAF0.02.recode.vcf.geno" file to your local laptop. You should be able to find this file in your home directory on the server from when we ran ADMIXTURE.
* Boot up R, set your working directory, and install the OutFLANK package from github:

```R
# Set your working directory to where you downloaded your results files:
setwd("~/github/PBIO381_srkeller_labnotebook/data/SNP_data/")

list.files() # Do you see your downloaded files there? If not, double check to make sure you've set your working directory to the right spot

install.packages("devtools")
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
install_github("whitlock/OutFLANK")

# ...and load the library. We'll also load the vcfR and adegenet libraries that we'll need later in the tutorial
library(OutFLANK)
library(vcfR)
library(adegenet)
```

* Next we'll import our SNP data in .geno format, and then transpose the data matrix (switch the rows and columns) to get it in the right orientation. 
* Also bring in the ssw meta-data that contains the info on disease status. Sort the data by individual ID, and then set the MM individuals disease "Trajectory" to NA. 

```R
ssw.geno <- read.fwf("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno", width=rep(1,24))
ssw.geno.trans <- t(ssw.geno)

ssw_meta <- read.table("ssw_healthloc.txt", header=T)
ssw_meta <- ssw_meta[order(ssw_meta$Individual),]
ssw_meta$Trajectory[which(ssw_meta$Trajectory=="MM")] = NA
```

* Now we're ready to run OutFLANK. This consists of the following steps:
  * Calculate Fst for each locus, based on divergence among the 3 disease groups (ssw_meta$Trajectory)
  * Trim back loci in the upper and lower 5% of the empirical Fst distribution, and eliminate low frequency SNPs (ones with expected heterozygosity, Hmin, < 0.1; equates to a MAF of ~0.05)
  * Re-fit the expected neutral Fst distribution to the data based on the presence of 3 groups in the data (NumberofSamples=3), and identify large Fst outliers in the upper tail based on a 10% false discovery rate (qthreshold=0.1). 
  * Plot the results: The empirical distribution of Fst (yellow bars) and the predicted neutral distribution (black line)

```R
OF_SNPs <- MakeDiploidFSTMat(ssw.geno.trans, locusNames=seq(1, 5317, by=1), popNames=ssw_meta$Trajectory)

OF_out <- OutFLANK(FstDataFrame=OF_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=3, qthreshold=0.1)

OutFLANKResultsPlotter(OF_out, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, titletext="Scan for local selective sweeps among disease groups: HH, HS, and SS")
```

* Take a look at your plot. Are there outlier loci that exceed the upper tail of this distribution?
* Identify your outliers and list which SNP loci they are, based on the rank order in the .geno matrix:

```
outliers <- which(OF_out$results$OutlierFlag=="TRUE")
print(outliers)
```

* Note that these numbers just refer to the rank order of the SNPs in your data matrix (from 1 to 5317). To figure out which genes (transcript IDs) these SNPs actually belong to, we need our vcf file. Make sure your vcf file is in the working directory of your R session, then read it in:

```
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
```

* We can extract the info on transcript ID and SNP position with the 'getFIX' command, and then ask for just this info associated with our outliers identified above:

```
vcfann <- as.data.frame(getFIX(vcf1))
vcfann[outliers,]
```

* Cool…it looks like we've got several promising candidate genes that show stronger than expected allele frequency divergence among our disease groups, consistent with selection imposed by SSW disease. 

* Coming up in our next coding session, we'll annotate our entire transcriptome using BLAST to public databases so we can assess what types of genes these are and what functions they have. 

* For now, we can do this manually by BLAST'ing a couple genes that we've identified as outliers:

  * Log back into the server and download the file that contains the full transcriptome reference sequences based on our RNASeq data:

    ```bash
    /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds
    ```

  * Open this file in your favorite light-weight text editor and search for one of the Trinity transcript ID's discovered by OutFLANK

  * Copy the resulting nucleotide sequences, and go to the NCBI BLAST homepage: https://blast.ncbi.nlm.nih.gov/Blast.cgi

  * Since these sequences are protein-coding transcripts, let's use the BLASTn algorithm, which translates the nucleotide sequence into amino-acids, and searches the NCBI nr database for protein matches

  * What did you find? What type of protein and its function are among the best hits?

  * It might be *really* interesting to do this for the outlier SNPs that you identified among disease groups using DAPC. Think back to the loading plots we made in the PopGenomics4 tutorial:

  ```R
  # Use the "droplevels" command to eliminate the last 2 individuals with "MM" and ensure there are only 3 disease categories; call this a new dataframe 'ssw_meta2'
  ssw_meta2 <- droplevels(ssw_meta[1:22,])

  # Now drop the last 2 individuals that are MM from the genlight SNP data. Note, the vcf object has an extra column in front, so you have to go to column 23 to get the first 22 individuals.
  gl1 <- vcfR2genlight(vcf1[,1:23])
  gl1$other <- as.list(ssw_meta2)

  # Run the DAPC using disease status to group samples
  disease.dapc <- dapc(gl1, pop=gl1$other$Trajectory, n.pca=7, n.da=3,
       var.loadings=T, pca.info=T, parallel=F)

  # Scatterplot of results
  scatter.dapc(disease.dapc, grp=gl1$other$Trajectory, legend=T)

  # Plot loci that contribute the most to distinguishing Healthy vs. Sick individuals (upper 0.1% of loadings)
  loadingplot(disease.dapc$var.load, 
              lab.jitter=1, 
              threshold=quantile(disease.dapc$var.load, probs=0.999))

  # What are the Trinity transcript ID's for these DAPC outliers?
  gl1$chromosome[which(disease.dapc$var.load>quantile(disease.dapc$var.load, 0.999))]
  ```

* Do any of the DAPC outliers overlap with the selected SNPs identified by OutFLANK? 

* What type of genes do you get for DAPC outliers when you use BLASTn?