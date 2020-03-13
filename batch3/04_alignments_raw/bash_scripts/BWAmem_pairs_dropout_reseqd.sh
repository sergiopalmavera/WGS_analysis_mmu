#!/bin/bash

REF=../../../reference_genome_ensembl
FASTQ=../../02_quality_trimming_adapter_removal/output
R1=I34772-L1_S19_L004_R1_001.corrected.fastq
R2=I34772-L1_S19_L004_R2_001.corrected.fastq
FLNM=I34772-L1_S19_L004.bam
time bwa mem -t 20 -M $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz $FASTQ/$R1 $FASTQ/$R2 | samtools view -@ 20 -bS - > ../output/$FLNM.bam 
