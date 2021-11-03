#!/bin/bash

pop=$1

fl_in=cohort_biallelicINDELs_VQSR99_PASS_withmissingnes.filtered.allrecords.${pop}.ALTfrq2 

fl_out=${fl_in/.ALTfrq2/.ALTfrq2_fixed_cols}

cat ../output/$fl_in | awk '{print $1 "\t" $2 "\t" $3}' | grep -v "^CHROM" > ../output/$fl_out

# the order of the columns is chr pos alt_fq
