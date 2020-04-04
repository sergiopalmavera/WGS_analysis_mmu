#!/bin/bash

hostname
date

chr=$1
db=$2
out=../output

mkdir $out/chr$chr
mkdir $out/chr$chr/tmp

#GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0 # version used until HaplotypeCaller
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0 # keep on using same version as for ConsolidateGVCFs (less bugs)
REF=../../reference_genome_ensembl

nm=$(basename $db)

$GATK/gatk GenotypeGVCFs \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V gendb://$db \
	-O $out/chr$chr/${nm}.vcf.gz \
	--tmp-dir=$out/chr$chr/tmp \
	--allow-old-rms-mapping-quality-annotation-data
rm -r $out/chr$chr/tmp
