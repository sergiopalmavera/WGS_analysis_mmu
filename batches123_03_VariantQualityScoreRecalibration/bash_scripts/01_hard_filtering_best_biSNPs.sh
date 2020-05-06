#!/bin/bash

# Each filter expression reads as i.e. "fail if QD < 20"

in_dir=../../batches123_02_GenotypeGVCFs_by_chr/output
in_vcf=cohort_biallelicSNPs.vcf.gz

out_dir=../output
out_vcf=cohort_biallelicSNPs_best.vcf

#GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

$GATK/gatk VariantFiltration \
	-R ../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $in_dir/$in_vcf \
	-filter-expression "MQ < 55.0" --filter-name "filter_MQ" \
	-filter-expression "FS > 0.1" --filter-name "filter_FS" \
	-filter-expression "QD < 25.0" --filter-name "filter_QD" \
	-filter-expression "MQRankSum < -0.5" --filter-name "filter_MQRankSumL" \
	-filter-expression "MQRankSum > 0.5" --filter-name "filter_MQRankSumR" \
	-filter-expression "ReadPosRankSum < -0.5" --filter-name "filter_ReadPosRankSumL" \
	-filter-expression "ReadPosRankSum > 0.5" --filter-name "filter_ReadPosRankSumR" \
	-filter-expression "SOR > 1.0" --filter-name "filter_SOR" \
	-O $out_dir/$out_vcf
