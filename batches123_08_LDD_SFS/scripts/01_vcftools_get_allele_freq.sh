#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

ls -1 ../../batches123_04_final_VCF/output/*.filtered.allrecords.*.vcf | for in_fl in $(cat)
do
	out_fl=$(basename $in_fl | sed 's/.vcf//')
	$vcftools/vcftools --freq --vcf $in_fl --out ../output/$out_fl
	wc -l ../output/$out_fl #check number of records
done
