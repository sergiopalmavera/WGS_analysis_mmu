#!/bin/bash

# Print hostname and data
hostname
date

chr=$1
chr_dir=../../batches123_01_ConsolidateGVFs_by_intervals/output/chr$chr

mkdir ./chr${1}_logs
logs=./chr${1}_logs

echo "# Joint Genotype consolidated gvcfs chromosome $1"

echo "## Total number of data stores:"
ls -1 $chr_dir | wc -l

for db in $(ls -1 $chr_dir)
do
	echo "## Submitting job for $db"
	nohup ./02_GenotypeGVCFs_template.sh $chr $chr_dir/$db &> ./$logs/GenotypeGVFs_${db}.out &
done

echo "### Process completed: all jobs submitted"

