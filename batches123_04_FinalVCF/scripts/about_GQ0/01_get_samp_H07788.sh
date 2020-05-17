#!/bin/bash


GATK4=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
GATK3=/home/fb4/palma-vera/FBN_HOME/Tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
REF=../../../reference_genome_ensembl/
VCF=cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.vcf

#java -jar $GATK3/GenomeAnalysisTK.jar -T SelectVariants \
#	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
#	-V ../../output/$VCF \
#	-o ./${VCF/.vcf/_H07788.vcf} \
#	-sn H07788-L1

$GATK4/gatk VariantsToTable \
	-V  ./${VCF/.vcf/_H07788.vcf} \
	-O ./${VCF/.vcf/_H07788.table} \
	-F CHROM -F POS -F REF -F ALT -GF GT -GF GQ

