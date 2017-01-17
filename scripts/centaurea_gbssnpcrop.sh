#!/bin/bash

# First script parses the raw reads contained in the fastq files. Needs the residual restriction site sequence after cut, which for ApeKI (used here for Centaurea libraries) is 'CWGC'

# Note, there is an issue with script 1 when the retriction recognition sequence contains an ambiguity code (see google groups discussion thread, "Wobble base in restriction site". The solution here was an edit to the *-1.pl script, where I edited out lines 25 and 26 (so no longer require to specify the recognition site in the comman line flag -enz1 or enz2) and also edited line 173 to manually provide the ApeKI recognition sequence, C[AT]GC. I saved this as a new script, GBS-SBO-CROP-1_ApeKI_edit.pl  

# Step 1: Screen raw fastq files, put all reads from multiple fastq into a single concatenated file, and change headers to barcode seq for each sample
# -------
# perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-1_ApeKI_edit.pl -d SE -b barcodes.txt -fq srkeller_GBS_20160812_C1_keller -s 1 -e 1 -enz1 CWGC

# Step 2: trim reads based on quality; remove adapter sequences
# -------
# perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-2.pl \
#	-d SE \
#	-fq srkeller_GBS_20160812_C1_keller \
#	-t 20 \
#	-ph 33 \
#	-ad /popgen/Trimmomatic-0.35/adapters/TruSeq3-SE.fa:2:30:10 \
#	-sl 4:30 \
#	-tr 30 \
#	-m 32 \
#	-l 30
 
# Step 3: Demultiplex reads into seaprate fastq files for each sample
# -------
# perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-3.pl \
#	-d SE \
#	-b /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/barcodes.txt \
#	-fq srkeller_GBS_20160812_C1_keller

# Before Step 4, need to install usearch into the system path (I did this into '/usr/bin'); don't need PEAR if only using single end reads

# Next error was thrown by useearch, which did not like the -fastaout option when building the MockRefGenome. After looking into the usearch documentation, I changed this flag in thte *4.pl script to -output and this solved the problem.

# Ahh, but another problem arose...usearch now runs out of memory when attempting to cluster the population centroids. This occurs even when randomly sampling the data, which is supposed to avoid the out of memory problem.

# After looking on the forum, I noticed this is an issue others have as well, and the open-source and 64-bit tool "vsearch" was offered as an alternative. I downloaded vsearch, installed to /usr/bin, and also downloaded an updated gbs-snp-crop perl script for step 3 that changes all clusteirng calls from "usearch" to "vsearch" under system calls during clustering.

# This would run, but would throw errors similar to above when trying to write output, due to non-recognition of the output flag, '-fastaout'. I editing these to now read '-output' based on the vsearch documentation. This new script is GBS-SNP-CROP-4_vsearch_edit.pl and works properly.

# Step 4: cluster reads and build the mock reference genome
# -------
#perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-4_vsearch_edit.pl \
#	-d SE \
#	-b /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/barcodes_nofail.txt \
#	-rl 100 \
#	-pl 32 \
#	-p 0.01 \
#	-id 0.93 \
#	-t 20 \
#	-MR Cent_C1_MockRef

# Clustered Mock Reference genome was made from randomly subsampling 6.4% of the total reads available, resulting in the file "UsearchIN.fasta" which contains 4.96 M reads.
# Total population clustering resulted in 1.75 M clusters with size (# reads) min=1 and max=124 (most are 1 or small).
# Clustered mock reference genome output as "Cent_C1_MockRef.MockRef_Genome.fasta", and should be used for read mapping in step 5


# Step 5: use bwa to map reads to mock reference; use samtools to pileup reads to call SNPs
# -------

# In the next step, we will use bwa to align the reads to the mock reference and samtools to call SNP positions for each sample

# perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-5.pl \
#	-d SE \
#	-b /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/barcodes_nofail.txt \
#	-ref Cent_C1_MockRef.MockRef_Genome.fasta \
#	-Q 30 \
#	-q 0 \
#	-f 0 \
#	-F 2308 \
#	-t 20 \
# 	-Opt "-m1 -s156"

# Step 6: parse mpileup output and produce SNP master matrix for downstream calling of SNP genotypes
# ---------

# perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-6.pl \
#	-b /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/barcodes_nofail.txt \
#	-out Centaurea_C1_SNP_master_matrix.txt

# Step 7: Filter SNPs and call genotypes

#perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-7.pl \
#	-in Centaurea_C1_SNP_master_matrix.txt \
#	-out Centaurea_C1_SNP_genotyping_matrix_50pct.txt \
#	-mnHoDepth0 5 \
#	-mnHoDepth1 20 \
#	-mnHetDepth 3 \
#	-altStrength 0.9 \
#	-mnAlleleRatio 0.1 \
#	-mnCall 0.5 \
#	-mnAvgDepth 4 \
#	-mxAvgDepth 200

# Step 8: file output in standard formats -- going to use the custom vcf script made available for this purpose:

perl /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/GBS-SNP-CROP-8_vcf_edit.pl \
	-in Centaurea_C1_SNP_genotyping_matrix_50pct.txt \
	-out Centaurea_C1_SNP_genotyping_matrix_50pct.vcf \
	-b /data/datashare/Centaurea/Centaurea_GBS/gbs_snp_crop/barcodes_nofail.txt \
	-formats vcf
	 



