## Update on installing R packakes needed for phyloseq microbiome analysis



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