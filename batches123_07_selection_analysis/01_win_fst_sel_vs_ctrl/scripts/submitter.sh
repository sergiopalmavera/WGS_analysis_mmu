#!/bin/bash

prefix=$1
pop1=$2
pop2=$3

script_nm=${prefix}_vcftools_fst_${pop1}_vs_${pop2}_per_SNP.sh

sed 's/pop1=/pop1='$pop1'/; s/pop2=/pop2='$pop2'/' template_fst_pop1_vs_pop2_per_SNP.sh > $script_nm

chmod +x $script_nm

nohup $script_nm &> ${script_nm/.sh/.out} &


