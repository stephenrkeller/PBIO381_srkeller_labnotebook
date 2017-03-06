#!/bin/bash

# Calling SNP variants using the reads2snp pipeline developed by N. Galtier and implemented in Gayral et al. 2014

##################################################
##  Step 1: Fix mate pairs in the SSW bam files ##
##################################################

cd /data/project_data/sam

for file in *.sam.bam
do
samtools fixmate "$file" >"./reads2snp/$file.fixmate.bam"
done

###################################
##  Step 2: Sort fixed bam files ##
###################################

cd reads2snp/

for file in *.fixmate.bam
do
sambamba_v0.6.0 sort -m 6G -t 8 -p --tmpdir=/data/project_data/sam/tmp2/ "$file" "$file.fixmate.sorted.sambada"
done

#############################
## Step 3: Remove PCR dups ##
#############################

for file in *.fixmate.sorted.sambada
do 
sambamba_v0.6.0 markdup "$file" "$file.fixmate.sorted.mrkdup.bam"
done

# This is now ready to pass to reads2snps for mapping

###################################################################
## Step 4: Create a list of the bam files for input to reads2snp ##
###################################################################

# The bam list file must be a two-colon file including one line per analyzed individual, such as:

# file1.bam   indiv1_name
# file2.bam   indiv2_name

# To make this list from the reads2snp/ directory where the sorted bam files are located:

ll *.bam >bamlist.txt
grep -o -P '\d\d\_.*bam' bamlist.txt >bamfiles.txt
grep -o -P '\d\d\_\d\-\d\d\_\w\_\d' bamfiles.txt >bamnames.txt
paste baminds.txt baminds2.txt >SSW_bamlist.txt

#############################################################
## Step 5: Call SNPs from bam files and reference assembly ##
#############################################################

/data/popgen/reads2snp_2.0/reads2snp_2.0.64.bin \
   -bam_list_file SSW_bamlist.txt  \
   -bam_ref_file /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds

