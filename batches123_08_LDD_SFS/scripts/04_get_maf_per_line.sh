#!/bin/bash

plink=/home/fb4/palma-vera/FBN_HOME/Tools/plink_v2

ls -1 ../../batches123_04_final_VCF/output/*.filtered.*.vcf | for in_fl in $(cat)
do
	echo "# Processing $in_fl"
	out_fl=$(basename $in_fl | sed 's/.vcf//')
	$plink/plink --vcf $in_fl --make-bed --out ../output/${out_fl}.plink
	$plink/plink --bfile ../output/${out_fl}.plink --freq --out ../output/${out_fl}.plink --allow-extra-chr
	printf "\n\n"
done



