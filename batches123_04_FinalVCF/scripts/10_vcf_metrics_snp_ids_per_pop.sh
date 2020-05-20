#!/bin/bash

pop=$1

REF=../../reference_genome_ensembl
PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

in_vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.vcf
pop_vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.${pop}.vcf
samps=../../sample_info/vcf_samples_$pop

echo "# subset main vcf for population (might contain fixed ref alleles)"
$bcftools/bcftools view --samples-file $samps $in_vcf -o $pop_vcf
printf "\n\n"

echo "# Making index"
$GATK/gatk IndexFeatureFile -F $pop_vcf 
printf "\n\n"

echo "# get metrics for pop"
java -jar $PICARD/picard.jar CollectVariantCallingMetrics \
	INPUT=$pop_vcf \
	OUTPUT=../metrics/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.${pop} \
	DBSNP=$REF/mus_musculus.vcf
printf "\n\n"

echo "# get SNPs for population"
$bcftools/bcftools view --samples-file $samps -i 'COUNT(GT="AA")>=1 || COUNT(GT="het")>=1' $pop_vcf -Ou | $bcftools/bcftools -f '%CHROM  %POS\n' > ${pop_vcf/.vcf/.SNPids}
printf "\n\n"


