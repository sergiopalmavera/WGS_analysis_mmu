#!/bin/bash

chr=5

# Path to gatk (absolute path - modify accordingly)
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0

# Path to batch1 gvcfs (from fastq obtained in 2018 with raw-avg-cvg >= 20) (Absolute paths - modify accodingly)
path_gvcf_batch1=/projekte/I2-SOS-FERT/05_HaplotypeCaller_GVCF/results

# Path to batch2 gvcfs (from fastqs obtained in 2018 with raw-avg-cvg < 20. They were reseqed to boost them to 20x) (Absolute paths - modify accodingly)
path_gvcf_batch2=/projekte/I2-SOS-FERT/05_HaplotypeCaller_GVCF2/results

# Path to batch3 gvfs (from fastqs obtained in 2020) (relative path)
path_gvcf_batch3=../../batch3/07_HaplotypeCaller/output

# Since batch2 samples were reseqed (boosted) they must be excluded from batch1
batch1_gvcfs=$(ls -1 $path_gvcf_batch1/*.gz | grep -v 'H07739\|H07753\|H07759\|H07774\|H07775\|H07777\|H07778\|H07784\|H07785\|H07787')

# Define batch2 and batch3 gvcfs 
batch2_gvcfs=$(ls -1 $path_gvcf_batch2/*.gz)
batch3_gvcfs=$(ls -1 $path_gvcf_batch3/*.gz | grep -v "I34772-L1_S63_L003.sorted.RG.dedup.bqsr.g.vcf.gz") #exclude "drop-out" sample (I34772)

# Combine gvcfs into one list
batch123_gvcfs=$(echo $batch1_gvcfs $batch2_gvcfs $batch3_gvcfs | for i in $(cat); do echo "-V $i"; done)

echo "# Number of gvcfs being consolidated"
for i in $batch123_gvcfs; do echo $i; done | grep gz | wc -l

$GATK/gatk GenomicsDBImport \
	$batch123_gvcfs \
	--genomicsdb-workspace-path ../chr$chr \
	--intervals $chr
