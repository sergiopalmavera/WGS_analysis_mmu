#/bin/bash

tab_snps=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ann.tab
tab_indels=cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.ann.tab

echo "# Extracting high or moderate effects from snp data"
time awk '{ if($6 == "HIGH" || $6 == "MODERATE") {print} }' ../output/$tab_snps > ../output/${tab_snps/.tab/.HIGH.MODERATE.tab}
printf "\n\n"

echo "# Extracting high or moderate effects from indel data"
time awk '{ if($6 == "HIGH" || $6 == "MODERATE") {print} }' ../output/$tab_indels > ../output/${tab_indels/.tab/.HIGH.MODERATE.tab}
printf "\n\n"
