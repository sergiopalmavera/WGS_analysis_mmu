#!/bin/bash
# Define date and print it
echo "Date of start"
date
DATE=`date +%Y%m%d_%H%M%S`

# Define path to PICARDTOOLS
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

# Creating dictionary from fasta file
java -jar $PICARD/picard.jar CreateSequenceDictionary \
	R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	O=$REF/Mus_musculus.GRCm38.dna.primary_assembly.dict

# Creating fasta index file
samtools faidx $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa
