#!/bin/bash

java -Xmx512M -jar /data/popgen/PGDSpider_2.0.9.0/PGDSpider2-cli.jar -inputfile ./SSW_all_biallelic.MAF0.02.recode.vcf -inputformat VCF -outputfile ./SSW_all_biallelic.MAF0.02.recode.vcf.geno -outputformat EIGENSOFT -spid ./vcf2admixture_SSW.spid
