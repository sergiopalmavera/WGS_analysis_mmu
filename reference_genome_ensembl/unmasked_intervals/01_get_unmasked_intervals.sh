#!/bin/bash

# based on: https://www.biostars.org/p/320184/

# Get intervals for regions not masked

UCSC_UTILS=/home/fb4/palma-vera/FBN_HOME/Tools/UCSC_utils/
BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2
SAMTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/samtools-1.9_installed/bin
FASTA=Mus_musculus.GRCm38.dna.primary_assembly.fa

$UCSC_UTILS/faToTwoBit ../$FASTA ${FASTA/.fa/.2bit}

$UCSC_UTILS/twoBitInfo -nBed ${FASTA/.fa/.2bit} N.bed

$SAMTOOLS/samtools faidx ../$FASTA

cut -f 1,2 ../${FASTA}.fai > ${FASTA/.fa/.genome}

$BEDTOOLS/bedtools.static.binary complement -i N.bed -g ${FASTA/.fa/.genome} > intervals_unmasked.bed 

