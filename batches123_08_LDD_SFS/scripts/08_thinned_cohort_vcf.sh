#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

vcf_path=../../batches123_04_FinalVCF/output

vcf_fl=cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.vcf

$vcftools/vcftools --vcf $vcf_path/$vcf_fl --recode --recode-INFO-all --thin 100000 --out ../output/${vcf_fl/.vcf/_thinned} 

echo "Number of lines in thinned vcf:"
grep -v '^#' ../output/${vcf_fl/.vcf/_thinned}* | wc -l

# Reference post: https://www.biostars.org/p/347796/	
# From the official manual page:
# "--thin <integer> Thin sites so that no two sites are within the specified distance from one another."

# separating snps by 100K bp results in ~23K SNPs in total. This is similar as what it was suggested by https://www.biostars.org/p/347796/
