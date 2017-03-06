#!/bin/bash

# Merge bam files from the same individual into 1 before snp calling

#cd /data/project_data/sam/

#ll *.bam | grep -o -P "\d\d\_" | sort -u >../snps/unique_inds.txt

cd /data/project_data/snps/reads2snps/

INDS=`cat unique_inds.txt`

for IND in $INDS

do

samtools merge $IND.merged.bam /data/project_data/sam/$IND*.bam

done

 
