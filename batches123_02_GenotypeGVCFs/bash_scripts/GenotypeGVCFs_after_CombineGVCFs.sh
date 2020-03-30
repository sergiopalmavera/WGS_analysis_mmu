#!/bin/bash

echo "# Printing host and date..."
hostname
date
printf "\n"

# Absoulte paths (modify accordingly)
#GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0 # version used until HaplotypeCaller
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0 # keep on using same version as for ConsolidateGVCFs (less bugs)

# Relative paths
REF=../../reference_genome_ensembl
OUT=../output
GVCF=/projekte/I2-SOS-FERT/GitHub/WGS_analysis_mmu/batches123_01_CombineGVFs/output

$GATK/gatk GenotypeGVCFs \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $GVCF/cohort.g.vcf \
	-O $OUT/cohort.vcf \
	--allow-old-rms-mapping-quality-annotation-data
