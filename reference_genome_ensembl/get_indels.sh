#!/bin/bash

# Extract indels from reference VCF

VCFTOOLS=~/FBN_HOME/Tools/vcftools_0.1.13/cpp/

$VCFTOOLS/vcftools --vcf ./mus_musculus.vcf --keep-only-indels --out ./mus_musculus_only_indels --recode --recode-INFO-all
