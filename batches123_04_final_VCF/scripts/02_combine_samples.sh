#!/bin/bash

# Based on this tutorial
# https://gatkforums.broadinstitute.org/gatk/discussion/53/combining-variants-from-different-files-into-one 
# https://gatkforums.broadinstitute.org/gatk/discussion/3328/using-selectvariants-to-select-for-multiple-expressions

# Note
# All VCF files are derived from the same VCF (VQSR95, filtered by SNP and PASS). Thus all vcf have the same sites and the union is the same as merging (inner joining)

GATK3=/home/fb4/palma-vera/FBN_HOME/Tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef #GenomeAnalysisTK.jar
out_vcf=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.vcf
in_vcfs=$(ls -1 ../TMP/*.vcf | for i in $(cat); do echo "-V $i"; done)

# CombineVariants is not yet added to GATK4 -- working fine!
java -Xmx1000g -jar $GATK3/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-R ../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	$(echo $in_vcfs) \
	-o ../output/$out_vcf

printf "\n"
echo "# N lines in last VCF"
grep -v "^#" ../output/$out_vcf | wc -l
printf "\n"
