#!/bin/bash

# Path to gatk (absolute path - modify accordingly)
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

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
batch3_gvcfs=$(ls -1 $path_gvcf_batch3/*.gz | grep -v I34772) #exclude "drop-out" sample (I34772)

# Combine gvcfs into one list
batch123_gvcfs=$(echo $batch1_gvcfs $batch2_gvcfs $batch3_gvcfs | for i in $(cat); do echo "--variant $i"; done)

echo "# Number of gvcf files to combine:"
echo $batch123_gvcfs | sed 's/--variant//g' | wc -w
printf "\n"

echo "# Running CombineGVCFs tool"
$GATK/gatk CombineGVCFs \
	-R ../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	$(echo $batch123_gvcfs) \
	-O ../output/cohort.g.vcf
printf "\n"

echo "# CombineGVCFs completed - check log file to corroborate"

