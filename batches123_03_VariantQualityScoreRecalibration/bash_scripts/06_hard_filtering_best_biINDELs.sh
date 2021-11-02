#!/bin/bash

in_dir=../../batches123_02_GenotypeGVCFs_by_chr/output
in_vcf=cohort_biallelicINDELs.vcf.gz

out_dir=../output
out_vcf=cohort_biallelicINDELs_best.vcf

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0

$GATK/gatk SelectVariants \
	-R ../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $in_dir/$in_vcf \
	-L ../data/truth_training_INDELs.intervals \
	-O $out_dir/$out_vcf
