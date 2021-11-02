#!/bin/bash

# Define paths
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl
OUT=/projekte/I2-SOS-FERT/05_HaplotypeCaller_GVCF2/results
SAMP=H07784

# Switch to folder with input files
cd /projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR2/BQSRbam

echo "# HaplotypeCaller (GVCF mode) for sample $SAMP"

# Define BAM-BQSR file
BAM=$(ls -1 *.bam | grep $SAMP)

echo "## BAM-BQSR file: $BAM" 

time $GATK/gatk HaplotypeCaller \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
 	-I $BAM \
 	-O $OUT/$SAMP.g.vcf.gz \
 	-ERC GVCF
