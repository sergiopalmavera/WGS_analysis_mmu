#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../reference_genome_ensembl

out_vcf=../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.vcf

$GATK/gatk SelectVariants \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V ../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.vcf \
	-L ../output/keep_indels_DU6_12_REST_15.intervals \
	-O ${out_vcf/.vcf/.TMP.vcf}

# some sites could loose the ALT allele after adding missingness. Keep sites having at least one ALT allele count
$bcftools/bcftools view -i 'COUNT(GT="AA")>=1 || COUNT(GT="het")>=1' ${out_vcf/.vcf/.TMP.vcf} -o $out_vcf 

rm ${out_vcf/.vcf/.TMP.vcf}
