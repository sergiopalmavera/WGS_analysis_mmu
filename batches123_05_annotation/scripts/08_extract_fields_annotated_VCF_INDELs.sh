#!/bin/bash

SnpEff=/home/fb4/palma-vera/FBN_HOME/Tools/SnpEff/snpEff
vcfEffOnePerLine=~/FBN_HOME/Tools/SnpEff/snpEff/scripts
IN=cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.ann.vcf

cat ../output/$IN | $vcfEffOnePerLine/vcfEffOnePerLine.pl | java -Xmx4g -jar $SnpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" > ../output/${IN/.vcf/.tab}             




