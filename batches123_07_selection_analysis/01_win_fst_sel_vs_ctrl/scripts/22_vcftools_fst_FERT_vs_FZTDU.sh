#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp
vcf_dir=../../../batches123_04_final_VCF/output
vcf_nm=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf
lines_dir=../../../sample_info
out=../output/FERT_vs_FZTDU

cat $lines_dir/vcf_samples_DUK $lines_dir/vcf_samples_DUC > ./tmp

time $vcftools/vcftools \
	--vcf $vcf_dir/$vcf_nm \
	--weir-fst-pop ./tmp \
	--weir-fst-pop $lines_dir/vcf_samples_FZTDU \
	--fst-window-size 50000 \
	--fst-window-step 25000 \
	--out $out 

rm ./tmp
