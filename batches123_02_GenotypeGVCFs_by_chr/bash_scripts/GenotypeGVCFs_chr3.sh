#!/bin/bash

chr=3

#GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0

REF=../../reference_genome_ensembl/

input=../../batches123_01_ConsolidateGVFs_by_intervals/output/chr$chr/DB

mkdir ../output/chr$chr

echo "# Number of databases for chr$chr"
ls -1 $input | wc -l

for db in $(ls -1 $input)
do
	echo "# $db"
	$GATK/gatk --java-options "-Xmx100G" GenotypeGVCFs \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-V gendb://$input/$db \
		-O ../output/chr$chr/${db}.vcf.gz \
		--allow-old-rms-mapping-quality-annotation-data
	printf "\n"
done
