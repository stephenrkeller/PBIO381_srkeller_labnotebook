#!/bin/bash

# Run ADMIXTURE to determine the number of genetic clusters in the SNP data, 
# and the ancestry proportions of each individual

# Remember the utility of 'for loops'?

for K in 1 2 3 4 5

do 

admixture -C 0.000001 -j10 --cv ./SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno $K \
| tee log${K}.out

done

grep -h CV log*.out >chooseK.txt



