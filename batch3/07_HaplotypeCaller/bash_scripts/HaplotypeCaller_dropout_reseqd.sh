#!/bin/bash

echo "# Working on:"
hostname
printf "\n"

# Absolute path (modify accordingly) 
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Relative paths
REF=../../../reference_genome_ensembl
OUT=../output
IN=../../05_alignments_addRG_dedup_bqsr/output

# bam file
b=I34772-L1_S19_L004.sorted.RG.dedup.bqsr.bam

# Define output file name
out_nm=$(echo $b | sed 's/.bam/.g.vcf.gz/')

echo "## HaplotypeCaller on bam $b ..."
time $GATK/gatk HaplotypeCaller \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-I $IN/$b \
	-O $OUT/$out_nm \
	-ERC GVCF
printf "\n"

echo "## bam $b completed (HaplotypeCaller)"
printf "\n"
