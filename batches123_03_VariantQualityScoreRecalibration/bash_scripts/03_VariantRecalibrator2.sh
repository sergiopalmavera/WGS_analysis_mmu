#!/bin/bash

#GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../reference_genome_ensembl
MGP=../../resource_mgp_sanger_REL-1505-SNPs_Indels
VARS=../../batches123_02_GenotypeGVCFs_by_chr/output
FL_DIR=../output

echo "# VariantRecalibrator SNP mode"
$GATK/gatk VariantRecalibrator \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $VARS/cohort_biallelicSNPs.vcf.gz \
	--resource TRUTH,known=false,training=true,truth=true,prior=12.0:$FL_DIR/cohort_biallelicSNPs_best.vcf \
	--resource sanger,known=false,training=true,truth=false,prior=10.0:$MGP/mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf \
	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$REF/mus_musculus.vcf \
	-an MQ \
	-an QD \
	-an FS \
	-an MQRankSum \
	-an ReadPosRankSum \
	-an SOR \
	-mode SNP \
	--target-titv 2.0 \
	-tranche 100.0 \
	-tranche 99.9 \
	-tranche 99.0 \
	-tranche 95.0 \
	-tranche 90.0 \
	--output $FL_DIR/cohort_biallelicSNPs_2.recal \
	-tranches-file $FL_DIR/cohort_biallelicSNPs_2.tranches \
	--rscript-file $FL_DIR/cohort_biallelicSNPs_plots_2.R

# TiTv modified from 2.15 to 2.0 according to MGPv5 (TiTv of known SNPs PASS, by applying CollectVariantCallingMetrics)
