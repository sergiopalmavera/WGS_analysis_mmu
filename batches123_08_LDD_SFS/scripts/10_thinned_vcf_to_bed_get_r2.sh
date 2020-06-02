#!/bin/bash

plink=/home/fb4/palma-vera/FBN_HOME/Tools/plink_v2

ls -1 ../output/*thinned.recode.*.vcf | for in_fl in $(cat)
do
	echo "# Processing $in_fl"
	out_fl=$(basename $in_fl | sed 's/.vcf//')
	$plink/plink --vcf $in_fl --make-bed --out ../output/${out_fl}.plink
	$plink/plink --bfile ../output/${out_fl}.plink --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 5000 -out ../output/${out_fl}.plink --allow-extra-chr
	printf "\n\n"
done

# "--ld-window-r2 0" ... output all r2 scores (min = 0)
# "--ld-window 999999" ... max number of snps in window 
# "--ld-window-kb 5000" ... window size: set to 5Mb (as suggested by Qanbari)
