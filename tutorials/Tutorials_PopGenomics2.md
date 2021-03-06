# P/BIO381 Tutorials

## Population Genomics 2: Estimating the diversity present within populations

### March 08, 2017

Where last we left off…...

We were filtering SNPs and looking at how SNPs may show signs of deviation from Hardy-Weinberg equilibrium, either as a heterozygote excess OR a deficit. Let's check this again, but now we're working with ALL the sequence libraries merged per individual. Yeah! Also, we're now going to use a compressed (gzipped) file to save some space. 

Let's first re-apply our filters, and zip up the resulting filtered output file. Then we can take a look at Hardy-Weinberg Equilibrium (HWE).

**PATH TO THE DATA:** 

```
/data/_project_data/snps/reads2snps/SSW_byind.txt.vcf.gz
```

*VCFtools filtering strategy (same as last session):*

```bash
$ vcftools --gzvcf SSW_byind.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.8  
$ gzip SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf
$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf --hardy
```

You can then bring the HWE output file (called "out.hwe") into **R** to take a look at which sites show deviation of observed from expected heterozygosity. Let's do this all together.

------------------------------

## Getting summary stats for downstream analysis and plotting in R##

In the last session, we got familiar with working with SNP data in VCF files, and doing some basic filtering. Now that we have a filtered SNP dataset that has high-quality sites in it, let's look at some different measures of genetic diversity in our sea star population. 

Keep in mind that the diversity of a population primarily reflects its **effective population size (Ne)**.  Ne is shaped by many different aspects of a species' life history and ecology (sex ratio, generation time, mating system, offspring number, and many more!) as well as the population's history (bottlenecks, population growth). As a result, looking at the diversity within populations (and comparing to other populations or species) is a critical step in understanding how ecology shapes genomes. 

Now: Let's take a more in-depth look at the diversity hidden within these sea star data. There are many different ways to look at the diversity within populations using SNPs. Here are some that we'll think about for today:


- **Nucleotide diversity (pi)**: The average number of pairwise differences between all individuals in the population. This is equivalent to the expected heterozygosity (2pq) for a bi-allelic locus.
- **Allele frequencies (*p* and *q*)**: What is the frequency of a given SNP? Usually defined in terms of the Major (common) and minor (rare) allele at each SNP locus.
- **Site Frequency Spectrum (SFS):** Also known as the "Allele Frequncy Spectrum". It is the histogram of allele frequencies across all SNP loci. In other words, how many loci are rare (say, frequency in the population between 0-0.1)? How many loci are common (0.4-0.5)? It turns out the shape of this distribution has an incredible amount of information in it…both about the population's demographic history (Ne, size changes) and also selection (purifying selection, positive selection)




We've been talking about **pi** in the last 2 papers, and seen that this relates to effective population size (Ne). Let's calculate **pi** first:

```bash
$ vcftools --gzvcf filenamevcf.gz --site-pi --out SSW_pi
```



### Diversity metrics based on subsetting your VCF files: ###

Many times we'll want to subset the total SNP dataset to analyze diversity in different groups. Say, compare allele frequencies in all the Healthy vs. Sick seastars. Or Intertidal vs. Subtidal. This is easy, you just need to create a separate text file containing which samples below to which groups so you can tell VCFtools how to split things up.

As an example, **let's compare the SNP frequencies for all loci between Healthy and Sick animals**. Perhaps there are some loci that contribute to a difference in pathogen susceptibility, which could be identified this way. Let's take a look.

First, you need to create text files containing the individual ID's for just Sick individuals. The following file has all the sample ID's in it that are part of your VCF file (SNPs=Y,N), along with Health Trajectory (HH,HS,SS,MM), and Location (INT,SUB): 

```
/data/project_data/snps/reads2snps/ssw_healthloc.txt
```

Use this file to get *just the **healthy** individual sample IDs*. 

How can we create these files in a clever, unix-y way?   HINT: grep!   ;)

Name it:

* **"H_OneSampPerInd.txt"**

Create another for *just the **sick** individuals*:

* **"S_OneSampPerInd.txt"**


We'll also want to remove all but the first column of data — the sample IDs. Here's a trick:

```
$ cut -f1 H_OneSampPerInd.txt >H_SampleIDs.txt
```

Do the same for Sick individuals.

Now that we have our individuals separated by disease status, we can call VCFtools to calculate allele frequencies separately for each population. This will require 2 separate calls to VCFtools.

**Allele Frequencies between Healthy and Sick individuals:**

```bash
$ vcftools --vcf filename.vcf --freq2 --keep H_SampleIDs.txt --out H_AlleleFreqs
```

```bash
$ vcftools --vcf filename.vcf --freq2 --keep S_SampleIDs.txt --out S_AlleleFreqs
```

Before we bring this into R for plotting, let's gather one more comparison between our groups. 



Recall that Wright's Fst measures allele frequency variance between groups, but standardizes the estimate based on the mean frequencies within groups. This gives a complimentary view to just comparing the raw frequencies.

**Fst between Healthy and Sick individuals:**

```bash
$ vcftools --vcf ~/filename.vcf --weir-fst-pop ~/H_SampleIDs.txt --weir-fst-pop ~/S_SampleIDs.txt --out HvS_OneSampPerInd
```

Now, we can import these datasets into R and make some plots to examine how the diversity varies in our dataset. We can get the data into R in one of 2 ways:

1. Download the results files to our laptops using scp or Fetch [MacOS] or Winscp [PC]
2. Stay on the server, and use the command-line version of R. This latter option can be more efficient if we want to do a quick look and make some simple plots, but isn't good for more complicated tasks. Here's how you would do R at the command-line on the server:

```R
[srkeller@pbio381 reads2snps]$ R

R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-redhat-linux-gnu (64-bit)

> fst <- read.table("HvS_OneSampPerInd.weir.fst",header=T)
> str(fst)
'data.frame':	442 obs. of  3 variables:
 $ CHROM : Factor w/ 111 levels "TRINITY_DN35598_c0_g1_TRINITY_DN35598_c0_g1_i1_g.5802_m.5802",..: 65 65 100 100 100 100 100 100 88 88 ...
 $ POS : int  4566 4665 978 1404 1722 3426 3729 3912 115 141 ...
 $ WEIR_AND_COCKERHAM_FST: num  0.0305 0.0085 0.0305 -0.0188 0.0732 ...
```

OK, looks like the data read in to R OK. Let's make a quick histogram of Fst and save as a pdf in our home directory on the server (~/)

```R
> pdf("HvS_Fst.pdf")
> hist(fst$WEIR_AND_COCKERHAM_FST, breaks=20, col="red")
> dev.off()
null device 
          1
```

* The first line creates the file and calls it by a name of our choosing
* The second line asks for a histogram plot of the data, specifying the x-axis to be broken into 20 bins and to plot each bar in red
* The last line turns the plotting device (dev) off. This tells R to save and close the plot. If you look in your home directoty on the server, you'll see you pdf waiting for you….



Let's also bring in the Allele Frequency results we got from VCFtools. Two interesting results to look for are:

1.   Calculating allele frequency differences between Healthy vs. Sick individuals for each SNP

	2. Calculate the Site Frequency Spectrum (SFS), which is simply a histogram of the allele frequencies across loci

    ​

When you're done, end your **R** session and return to the command line:

```bash
> quit()
Save workspace image? [y/n/c]: n
[srkeller@pbio381 ~]$ 
```



------------------------------------------

###Comparing Sea Star nucleotide diversity and piN/piS to the sample of metazoans that Romiguier et al. (2014) report.###

Romiguier et al. report on some very intriguing associations between species life history traits, nucleotide diversity at synonymous sites (piS), and the ratio of piN/piS (where piN is nucleotide diversity at nonsynonymous site). 

![Romiguier_Figure2](http://www.nature.com/nature/journal/v515/n7526/images/nature13685-f2.jpg)



What do our sea star data have to say about this? Or more importantly: ***where do sea stars fall on the genomic diversity ~ life history continuum?***



Estimating piS and piN on the entire dataset will take some time. Let's try and set this up at the end of the day and let it run. We'll use the piNpiS program from Gayral et al. (2013) to run this. We only need a single input file, which is a FASTA formatted sequence file that is output from **reads2snps**

```bash
$ cd /data/project_data/snps/reads2snps
[srkeller@pbio381 reads2snps]$ /data/popgen/dNdSpiNpiS_1.0 -alignment_file=SSW_byind.txt.fas -ingroup=sp
```

While we wait for that to chug along, we can look at the output from the smaller VCF file run on just 1 sample library per individual. The output gives several files; here are the important ones:

* Detailed results file with nucleotide diversity calculated for each gene:

  * ```
    /data/project_data/snps/reads2snps/SSW_bamlist.txt.out
    ```

* Overall summary of results, including mean (and confidence intervals) of diversity when averaged across all expressed genes:

  * ```bash
    /data/project_data/snps/reads2snps/SSW_bamlist.txt.sum
    ```



To compare the mean values across genes to Romiguier's data, we need to get their estimates and combine them with our estimates. I've downloaded Table S3 to our server and saved as a common separated (.csv) file:

```bash
/data/project_data/snps/reads2snps/Romiguier_nature13685-s3.csv
```



