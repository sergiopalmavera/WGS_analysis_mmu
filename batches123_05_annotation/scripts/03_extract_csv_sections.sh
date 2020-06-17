#!/bin/bash

#had to do this by hand!

fl_nm=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.summary.csv

cat ../output/$fl_nm | sed -n '4,18p' > ../output/${fl_nm/.csv/_summary_table.csv}
cat ../output/$fl_nm | sed -n '22,42p' > ../output/${fl_nm/.csv/_change_rate_by_chr.csv}
cat ../output/$fl_nm | sed -n '46,47p' > ../output/${fl_nm/.csv/_vars_by_type.csv}
cat ../output/$fl_nm | sed -n '51,55p' > ../output/${fl_nm/.csv/_eff_by_impact.csv}
cat ../output/$fl_nm | sed -n '59,62p' > ../output/${fl_nm/.csv/_eff_by_func_class.csv}
cat ../output/$fl_nm | sed -n '68,87p' > ../output/${fl_nm/.csv/_cts_by_eff.csv}
cat ../output/$fl_nm | sed -n '91,102p' > ../output/${fl_nm/.csv/_cts_by_geno_region.csv}
