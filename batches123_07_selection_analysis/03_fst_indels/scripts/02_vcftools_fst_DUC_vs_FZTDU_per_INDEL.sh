#!/bin/bash
pop1=DUC
pop2=FZTDU
vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp
vcf_dir=../../../batches123_04_final_VCF/output
vcf_nm=cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.vcf
lines_dir=../../../sample_info
out=../output/${pop1}_vs_${pop2}

time $vcftools/vcftools --vcf $vcf_dir/$vcf_nm --weir-fst-pop $lines_dir/vcf_samples_$pop1 --weir-fst-pop $lines_dir/vcf_samples_$pop2 --out $out 
