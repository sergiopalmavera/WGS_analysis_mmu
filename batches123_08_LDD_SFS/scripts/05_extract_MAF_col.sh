#!/bin/bash

for maf_fl in  ../output/*.plink.frq
do
	awk '{print $5}' $maf_fl > ${maf_fl/.frq/.maf}
done
