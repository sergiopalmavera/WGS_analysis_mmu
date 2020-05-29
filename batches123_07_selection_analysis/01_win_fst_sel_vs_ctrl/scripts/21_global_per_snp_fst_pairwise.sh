#!/bin/bash

ls -1 *_per_SNP.out | for i in $(cat)
do
	cont=$(echo $i | sed 's/vcftools_fst_//; s/_per_SNP.out//; s/[[:digit:]][[:digit:]]_//')
	global_fst=$(grep "Weir and Cockerham mean Fst estimate:" $i | cut -d':' -f2 | sed 's/ //')
	echo $cont $global_fst
done

