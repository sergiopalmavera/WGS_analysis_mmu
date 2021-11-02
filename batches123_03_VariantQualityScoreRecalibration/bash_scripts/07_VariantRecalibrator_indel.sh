#!/bin/bash

#https://github.com/broadgsa/gatk/blob/master/doc_archive/faqs/Which_training_sets___arguments_should_I_use_for_running_VQSR%3F.md
#"Note that indels use a different set of annotations than SNPs. Most annotations related to mapping quality have been removed since there is a conflation with the length of an indel in a read and the degradation in mapping quality that is assigned to the read by the aligner. This covariation is not necessarily indicative of being an error in the same way that it is for SNPs."

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../reference_genome_ensembl
MGP=../../resource_mgp_sanger_REL-1505-SNPs_Indels
VARS=../../batches123_02_GenotypeGVCFs_by_chr/output
FL_DIR=../output

$GATK/gatk VariantRecalibrator \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $VARS/cohort_biallelicINDELs.vcf.gz \
	--resource TRUTH,known=false,training=true,truth=true,prior=12.0:$FL_DIR/cohort_biallelicINDELs_best.vcf \
	--resource sanger,known=false,training=true,truth=false,prior=10.0:$MGP/mgp.v5.merged.indels.dbSNP142.normed.PASS.vcf.gz \
	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$REF/mus_musculus.vcf \
	-an QD \
	-an FS \
	-an ReadPosRankSum \
	-an SOR \
	-mode INDEL \
	--target-titv 2.0 \
	-tranche 100.0 \
	-tranche 99.9 \
	-tranche 99.0 \
	-tranche 95.0 \
	-tranche 90.0 \
	--output $FL_DIR/cohort_biallelicINDELs.recal \
	-tranches-file $FL_DIR/cohort_biallelicINDELs.tranches \
	--rscript-file $FL_DIR/cohort_biallelicINDELs_plots.R

# TiTv modified from 2.15 to 2.0 according to MGPv5 (TiTv of known SNPs PASS, by applying CollectVariantCallingMetrics)
