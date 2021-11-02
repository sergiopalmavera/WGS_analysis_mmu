#!/bin/bash

pop=$1
fl=../../sample_info/vcf_samples_$pop

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
vcftools=~/FBN_HOME/Tools/vcftools_0.1.13/cpp

in_vcf_path=../../batches123_03_VariantQualityScoreRecalibration/output
in_vcf=cohort_biallelicINDELs_VQSR99_PASS.vcf

for samp_nm in $(cat $fl); do
#samp_nm=I34713
	echo "# Sample $samp_nm"

	# Extract sample from main VCF
	$bcftools/bcftools view -s $samp_nm $in_vcf_path/$in_vcf -o ../TMP/${in_vcf/.vcf/_${samp_nm}.vcf}

	# Extract DP field 
	$bcftools/bcftools query -f '[%DP ]\n' ../TMP/${in_vcf/.vcf/_${samp_nm}.vcf} -o ../TMP/${samp_nm}.DP

	# Calculate mean DP
	mean_DP=$(cat ../TMP/${samp_nm}.DP | tr ' ' \\n | awk '{x+=$0; next} END{print x/NR}')
	echo "mean DP: $mean_DP"

	# Calculate sd DP (formula from https://stackoverflow.com/questions/18786073/compute-average-and-standard-deviation-with-awk)
	sd_DP=$(cat ../TMP/${samp_nm}.DP | tr ' ' \\n | awk '{sum+=$0;a[NR]=$0}END{for(i in a)y+=(a[i]-(sum/NR))^2;print sqrt(y/(NR-1))}')
	echo "Standard Deviation: $sd_DP"

	# Calculate maxDP (according to post: https://www.biostars.org/p/265782/)
	max_DP=$( echo "$mean_DP + 3*$sd_DP" | bc )
	echo "max DP: $max_DP"

	# Add missingness to sample VCF
	$vcftools/vcftools --vcf ../TMP/${in_vcf/.vcf/_${samp_nm}.vcf} --minGQ 20 --minDP 4 --maxDP $max_DP --recode --recode-INFO-all --out ../TMP/${in_vcf/.vcf/_${samp_nm}_missingness}

	# Remove files to free space
	rm ../TMP/${in_vcf/.vcf/_${samp_nm}.vcf} ../TMP/${samp_nm}.DP
done



