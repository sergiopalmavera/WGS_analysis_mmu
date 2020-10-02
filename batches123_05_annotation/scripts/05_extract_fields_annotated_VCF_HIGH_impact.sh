#/bin/bash

awk '{ if($6 == "HIGH") {print} }' ../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ann.tab > ../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ann.HIGH.tab
