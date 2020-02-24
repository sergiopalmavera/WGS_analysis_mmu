#!/bin/bash

echo "# Working on:"
hostname
printf "\n"

# Define samples in batch
FROM=44
TO=47

# Absolute path (modify accordingly) 
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Relative paths
REF=../../../reference_genome_ensembl
OUT=../output
IN=../../05_alignments_addRG_dedup_bqsr/output

# Define samples (reconstruct bam file names from raw fastq files)
bam_nms=$(ls -1 ../../01_quality_control/output/*.zip | for i in $(cat); do basename $i | sed 's/_fastqc.zip/.sorted.RG.dedup.bqsr.bam/;s/_R[1-2]_001//' ; done | sort | uniq)

# Define samples in batch
bam_files=$(for b in $bam_nms; do echo "$IN/$b" ;done | sed -n $FROM,${TO}p)

echo "# Processing files"
for b in $bam_files; do echo $b; done
printf "\n"

time for b in $bam_files;
do
	# Define output file name
	out_nm=$(basename $b | sed 's/.bam/.g.vcf.gz/')

	echo "## HaplotypeCaller on bam $b ..."
	time $GATK/gatk HaplotypeCaller \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $b \
		-O $OUT/$out_nm \
		-ERC GVCF
	printf "\n"

	echo "## bam $b completed (HaplotypeCaller)"
	printf "\n"
done

echo "# Done processing files:"
for b in $bam_files; do echo $b; done
