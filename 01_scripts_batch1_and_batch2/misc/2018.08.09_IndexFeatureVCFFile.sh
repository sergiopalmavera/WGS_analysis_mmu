#!/bin/bash
# Define date and print it
echo "Date of start"
date
DATE=`date +%Y%m%d_%H%M%S`

# Define path to GATK
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0


# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

$GATK/gatk IndexFeatureFile \
	-F $REF/mus_musculus.vcf
