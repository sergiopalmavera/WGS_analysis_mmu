#!/bin/bash

from_sample=$1
to_sample=$2

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
vcf_path=../../batches123_04_final_VCF/output
vcf_fl=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf
out_fl=samp${from_sample}_to_samp${to_sample}_with_rec_rate.roh

sed -n "$from_sample,${to_sample}p" ../../sample_info/vcf_samples.txt > ./TMP/samples_${from_sample}_${to_sample}.txt

echo "# Processing samples:"
cat ./TMP/samples_${from_sample}_${to_sample}.txt
printf "\n"
date
printf "\n"
time $bcftools/bcftools roh --rec-rate 0.51e-8 --AF-file ../output/${vcf_fl/.vcf/.AF.tab.gz} $vcf_path/$vcf_fl -o ../output/$out_fl -S ./TMP/samples_${from_sample}_${to_sample}.txt

# previously I ran roh in default mode without specifzing the recombination rate.
# i could not find the defaults, but it appears that it is 1e-8, based on the paper ('Mountain gorilla genomes reveal the impact of long-term population decline and inbreeding') that preceeded the publication of bcftools roh
# 'assumed constant here and equal to 1 x 10^-8 per base pair'

# I found that the recombination rate for mouse is:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8328415/ 
## 'For simplicity, we assumed a recombination rate of 0.51 cM/Mb, equal to the genomic average for house mice (Cox et al., 2009)'

# This corresponds to 0.51% chance of recombination per Mb
# --> 0.51%/Mb --> 0.51*10^-2 / 1*10^6 --> 0.51*10^-2*10^-6 --> 0.51*10^-8 --> 0.51e-8 
# This is half than the default.

# This is a good link to transform cM/Mb: https://www.biostars.org/p/285449/

# another reference about rec rate in mouse: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC383296/
# 'the genome-wide average recombination rates are 0.555 cM/Mb, 0.528 cM/Mb, and 1.20 cM/Mb for rat, mouse, and human, respectively. A more accurate measure including only distances measured between placed markers (i.e., not counting portions of chromosomes before and after the first and last markers), gives genome-wide estimates of 0.60 cM/Mb for rat, 0.56 cM/Mb for mouse, and 1.26 cM/Mb for human.'
