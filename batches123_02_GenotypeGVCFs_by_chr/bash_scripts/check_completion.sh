#!/bin/bash

script=$1
head -2 $script
printf "\n"

echo "# Number of data bases joint-genotyped"
grep "org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done" $script | wc -l
printf "\n"

echo "# Printing the last lines of output file"
grep -n "org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done" $script | cut -d':' -f1 | for i in $(cat);
do
	xx=$(( i - 3 ))
	ii=$(( i + 2 ))
	sed -n "$xx,${ii}p" $script
	printf "\n"
done

