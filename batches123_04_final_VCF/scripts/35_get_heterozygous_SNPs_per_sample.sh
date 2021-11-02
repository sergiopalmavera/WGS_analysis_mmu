#!/bin/bash

pop=$1
vcf=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf
bcftools=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin


while read s 
do
	echo "# Processing $s"
	out_fl=../output/${vcf/.vcf/.${s}.het.tab}

	time $bcftools/bcftools view -Ou -s $s ../output/$vcf | $bcftools/bcftools query -f '[%CHROM %POS %GT\n]' | grep -Ev '0/0|1/1|\./\.' > ../output/$out_fl

done < ../../sample_info/vcf_samples_${pop}
