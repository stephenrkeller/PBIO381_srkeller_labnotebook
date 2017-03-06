# P/BIO381 Tutorials

## Population Genomics 1: Intro to working with SNP data in variant call format (vcf), and manipulation with 'vcftools'

### March 06, 2017

When doing population genomics on large genome-wide or transcriptome-wide datasets, we generally want to work with files that contain just the polymoprhic sites and omit sites that are fixed. But there's also a lot of metadata from our assembly (remember those good 'ol sam files?) that will be important for analyzing the SNP data downstream. 

*These are things like:*

- **Position**: Where is the SNP located within a contig or chromosome of the reference assembly?
- **Alleles**: What are the alleles present at a given SNP? Are there only 2, or are there more? Are they single-nucleotide differences?
- **Depth**: How many reads cover a given SNP? How many reads were observed for each allele?
- **Genotype Quality (GQ):**  How confident are we that we're calling the correct genotype (ex., AA, AT, or TT)?
- **Sample Names:** Where are the data for each individual sample?



As usual, the community has converged on a common standard to represent these large and sometimes complex SNP data files. It is known as the Variant Call Format, or VCF. Here's a link to the description of what each field in a VCF file means:  [VCF version 4.3 file definition]([hts-specs](https://github.com/samtools/hts-specs)/**VCFv4.3.pdf**)



We'll be working with vcf files a lot as we conduct the population genomics section of the course. The first step in learning how to work with these files is to use a program called **VCFTools** for parsing your data file into just those samples and sites of interest, and to calculate diversity stats on these.



The manual page for VCFtools is an excellent resource! [The latest version is here.](https://vcftools.github.io/man_latest.html) 



## Basic Syntax and Usage ##

We're going to use VCFtools to examine the effects of different filtering strategies on the number of SNPs we get and their quality. The first step is seeing if VCFtools likes our file format, and getting some basic info on the # of SNPs and samples:

```bash
$ vcftools --vcf filename.vcf
```

This will return some basic info that should match of general expectations of sample size. It's always good to have a reality check on your file format.



Let's try filtering out positions that are likely to be errors in the sequencing or genotyping process. For now, let's just identify how many SNPs would pass each filter without actually changing the datafile at all. Then, we can decide what combination of filters we may want to implement.



### Record for each of the following steps the number of SNPs (aka sites) that would be make it through each filter:###



* *Biallelic vs. multi-allelic SNPs:*  Get rid of sites with >2 alleles. 
  * Rationale: When looking at diversity within species, it's very rare to have mutations occur at the same position. So, SNPs with >2 alleles probably reflect sequence or mapping errors.

```bash
$ vcftools --vcf filename.vcf --max-alleles 2
```



* *Read depth per SNP:* If an individual has <10 reads mapped at a given SNP,  convert to missing data
  * Rationale: To have confidence that the SNP alleles for an individual have been sampled, you need read depth (the more the better). This is the coin flipping analogy…how many times do you need to flix the coin to be confident that it is fair, with an equal chance of heads & tails.

```bash
$ vcftools --vcf filename.vcf --minDP 8
```



* *Genotype Quality:* Exclude all genotype calls below a quality threshold of Q20 (corresponds to 1 in 100 probability of being incorrect)
  * Rationale: How confident are you that your homozygous call isn't really a heterozygote and you just missed seing the 2nd allele? The GQ scores gives you a PHRED-scaled measure of confidence in the genotype call being correct, based on the multinomial probability of getting the 3 genotypic class, given the number of reads covering each allele.

```bash
$ vcftools --vcf filename.vcf --minGQ 20
```



- *Missing data across individuals:* Get rid of sites where  fewer than 80% of our samples have data. 
  - Rationale: Missing data is a problem for any analysis, and population genetic statistics can behave oddly (i.e.. become biased) when a lot of individuals are missing data for a given SNP. 

```bash
$ vcftools --vcf filename.vcf --max-missing 0.8
```



### Combining filters: ###

Now, it's time to combine filters instead of applying them one at a time. **NOTE:** VCFtools processes the filter requests in the order that you give it at the command-line. This is a key point, and means that if you apply the same filters in different orders, you will likely get different results!

* I recommend the following order:   biallelic filter>depth>GQ>missingness
* To output the resulting filtered data as a new vcf file, add the "—recode —out outfilename" to the end of the command. 

```bash
$ vcftools --vcf filename.vcf --max-alleles 2 --minDP 10 --minGQ 20 --max-missing 0.8 --recode --out outputfilename
```









​

##Let's review what we've learned so far…##

 

******************

