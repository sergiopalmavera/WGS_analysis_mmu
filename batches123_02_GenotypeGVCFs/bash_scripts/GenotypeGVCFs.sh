#!/bin/bash

# Absoulte paths (modify accordingly)
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Relative paths
GVCF=../../batches123_01_CombineGVFs/output
REF=../../reference_genome_ensembl
OUT=../output

$GATK/gatk GenotypeGVCFs \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $GVCF/cohort.g.vcf \
	-O $OUT/joint_genotype_cohort.vcf 
