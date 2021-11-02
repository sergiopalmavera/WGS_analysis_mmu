#!/bin/bash

script=$1
echo "# $script"
chr=$(grep "chr=" $script | sed 's/chr=//')
dir_logs=../output/chr$chr/logs
dir_intervals=../output/chr$chr/intervals

cond=$( echo $script | grep "_pt" | wc -l)
cond2=$( echo $script | grep "_from" | wc -l)

if [ $cond == 1 ]
then
	pt=$(echo $script | cut -d'_' -f5 | sed 's/pt// ; s/.sh//')
	intervals=$dir_intervals/chr${chr}_pt$pt.bed	
else
	if [ $cond2 == 1 ]
	then
		pt=$(echo $script | cut -d'_' -f5 | sed 's/.sh//')	
		intervals=$dir_intervals/chr${chr}_$pt.bed
	else
		intervals=$dir_intervals/chr${chr}.bed
	fi
fi

while read l
do
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}')
	echo "$dir_logs/ConsolidateGVCFs_chr${nm}.out ($(echo $l | awk '{print $3-$2}')bp)"
	grep "GenomicsDBImport done" $dir_logs/ConsolidateGVCFs_chr${nm}.out
	grep -E "ERROR|error|Error" $dir_logs/ConsolidateGVCFs_chr${nm}.out
done < $intervals > ./TMP/tmp

echo "## Number of jobs/intervals completed"
cat ./TMP/tmp | grep "GenomicsDBImport done" | wc -l

echo "## Number of jobs/intervals submitted"
xx=$(echo $script | sed 's/from// ; s/to/_/ ; s/.sh//')
xx1=$(echo $xx | cut -d'_' -f5)
xx2=$(echo $xx | cut -d'_' -f6)
echo $(( ($xx2 - $xx1)+1 ))

cat ./TMP/tmp
printf "\n"


