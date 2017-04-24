### Update on creating your bio file for testing functional enrichment of your microbiome using Picrust

Turns out that there are 2 different types of 'biom' formatted files — an older version called "JSON", and a newer version called "HDF5". By default, PICRUST outputs the new HDF5 version, but the libraries we will use in R will not recognize this! So, to make all the packages work together, we need to convert from out text formatted PICRUST output to JSON format, and then bring the JSON biom file into R.

We'll do this on the class server. Follow these steps:

1. cd to the directory with your *L3.biom file that you ouput as the last step of your PICRUST analysis; for a refresher, see [the tutorial](https://adnguyen.github.io/2017_Ecological_Genomics/Tutorial/2017-04-19_picrust.html)
2. while in this directory, execute the following command:

```
biom convert -i metagenome_predictions.L3.txt -o metagenome_predictions.L3.json.biom --table-type="OTU table" --to-json

```

It should only take a few seconds. Then, download the new *.json.biom file to your desktop and follow Melanie's R script at the bottom of the [tutorial](https://adnguyen.github.io/2017_Ecological_Genomics/Tutorial/2017-04-19_picrust.html) to test for functional enrichment of the micro biome. NOTE, if you need help getting the packages installed to run this analysis, see below update:



----------

## Update on installing R packages needed for phyloseq microbiome analysis



Install R package RJSONIO form within R Studio

```R
install.packages("RJSONIO")
```

Here's what you should see:

```
trying URL 'http://cran.rstudio.com/bin/macosx/mavericks/contrib/3.3/RJSONIO_1.3-0.tgz'
Content type 'application/x-gzip' length 1279768 bytes (1.2 MB)
==================================================
downloaded 1.2 MB
```

Now, we need to install the 'biom' package. This is no longer available through the normal CRAN repository, so we need to download it first manually from the source website. You can use the following link to download the file to your hardrive (make sure you remember *where* you've chosen to download it):

https://cran.r-project.org/src/contrib/Archive/biom/biom_0.3.12.tar.gz

Don't try and open this file. Instead, go back to R studio and type the following command:

```R
install.packages("~/Downloads/biom_0.3.11.tar.gz", repos = NULL, type = "source")
```

Here's what you should see:

```
* installing *source* package ‘biom’ ...
** package ‘biom’ successfully unpacked and MD5 sums checked
** R
** inst
** preparing package for lazy loading
Creating a generic function for ‘nrow’ from package ‘base’ in package ‘biom’
Creating a generic function for ‘ncol’ from package ‘base’ in package ‘biom’
Creating a generic function for ‘rownames’ from package ‘base’ in package ‘biom’
Creating a generic function for ‘colnames’ from package ‘base’ in package ‘biom’
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded
* DONE (biom)
```



You should now be able to type:

```
library(biom)
```

…and continue with the rest of the R phyloseq tutorial for our 16S data posted by Melanie.