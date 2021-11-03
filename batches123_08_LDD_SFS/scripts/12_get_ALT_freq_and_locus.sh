#!/bin/bash

fls=$(ls -1 ../output/*.frq | grep -v "plink")

for fl in $fls
do
	echo $fl
	awk '{print $6}' $fl | awk -F ":" '{print $2}' > ./TMP/$(basename $fl).altfrq.TMP
	awk '{print $1 "\t" $2}' $fl > ./TMP/$(basename $fl).locus.TMP
	paste --delimiters=' ' ./TMP/$(basename $fl).locus.TMP ./TMP/$(basename $fl).altfrq.TMP > ${fl/frq/ALTfrq2}
	printf "\n"
done
