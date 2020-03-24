#!/bin/bash

chr=

# Absoulte paths (modify accordingly)
#GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0 # version used until HaplotypeCaller
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0 # keep on using same version as for ConsolidateGVCFs (less bugs)

# Relative paths
DB=../../batches123_01_ConsolidateGVFs/chr$chr
REF=../../reference_genome_ensembl
OUT=../output
TMP_DIR=./TMP/chr$chr

mkdir $TMP_DIR

$GATK/gatk GenotypeGVCFs \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V gendb://$DB \
	-O $OUT/chr${chr}.vcf.gz \
	--tmp-dir=$TMP_DIR \
	--allow-old-rms-mapping-quality-annotation-data

rm -r $TMP_DIR 
