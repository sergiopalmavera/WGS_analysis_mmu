#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../reference_genome_ensembl
VARS=../../batches123_02_GenotypeGVCFs_by_chr/output
FL_DIR=../output

pct_sens=99.9
echo "# Sensitivity to truth set = $pct_sens"
printf "\n"

echo "# ApplyRecalibrator SNP mode"
$GATK/gatk ApplyVQSR \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $VARS/cohort_biallelicSNPs.vcf.gz \
	-O $FL_DIR/cohort_biallelicSNPs_VQSR${pct_sens}.vcf \
	-ts-filter-level $pct_sens \
	--tranches-file $FL_DIR/cohort_biallelicSNPs_2.tranches \
	--recal-file $FL_DIR/cohort_biallelicSNPs_2.recal \
	-mode SNP
