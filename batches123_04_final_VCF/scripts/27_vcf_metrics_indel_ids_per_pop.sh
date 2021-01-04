#!/bin/bash

pop=$1

REF=../../reference_genome_ensembl
PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

in_vcf=../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.vcf
pop_vcf=../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.${pop}.vcf
samps=../../sample_info/vcf_samples_$pop

echo "# subset main vcf for population (exclude ref alleles)"
$bcftools/bcftools view --samples-file $samps $in_vcf -Ou | $bcftools/bcftools view -i 'COUNT(GT="AA")>=1 || COUNT(GT="het")>=1' -o $pop_vcf
printf "\n\n"

echo "# Number of INDELs in pop:"
grep -v '^#' $pop_vcf | wc -l
printf "\n\n"

echo "# Making index"
$GATK/gatk IndexFeatureFile -F $pop_vcf 
printf "\n\n"

echo "# get metrics for pop"
java -jar $PICARD/picard.jar CollectVariantCallingMetrics \
	INPUT=$pop_vcf \
	OUTPUT=../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.${pop} \
	DBSNP=$REF/mus_musculus.vcf
printf "\n\n"

echo "# get INDELs for population"
$bcftools/bcftools query -f '%CHROM  %POS\n' $pop_vcf > ${pop_vcf/.vcf/.INDELids}
printf "\n\n"

echo "# Number of INDEL ids"
wc -l ${pop_vcf/.vcf/.INDELids}

