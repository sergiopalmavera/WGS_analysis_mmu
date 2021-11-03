#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

for vcf in ../../../batches123_04_final_VCF/output/*filtered.allrecords.*.vcf
do
	echo "# Processing $vcf"
	out_fl=$(basename $vcf)
	$vcftools/vcftools --vcf $vcf --window-pi 50000 --window-pi-step 25000 --out ../output/${out_fl/.vcf/}
	echo "## Done"
	printf "\n"
done
