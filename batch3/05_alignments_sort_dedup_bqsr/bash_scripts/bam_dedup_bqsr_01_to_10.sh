#!/bin/bash

# Define samples in batch
FROM=1
TO=10

# Absolute paths (modify accordingly)
PICARD=~/FBN_HOME/Tools/picard_2.18.11
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Relative paths
bam_files_path=../../04_alignments_raw/output
REF=../../../reference_genome_ensembl
TMP=../TMP

# Define samples
bam_files=$(ls -1 $bam_files_path/*.sorted.bam | sed -n $FROM,${TO}p)

echo "# Processing files"
for b in $bam_files; do echo $b; done
printf "\n"

time for b in $bam_files
do
	bam=$(basename $b)
	
	echo "## BAM file $bam"
	mkdir $TMP/${bam}.TMP
	
	echo "## Marking duplicates ..."
	time java -jar $PICARD/picard.jar MarkDuplicates \
		I=$b \
		O=$TMP/${bam}.TMP/${bam/.bam/.dedup.bam} \
		M=..output/${bam/.bam/.dedup.metrics.txt} \
		TMP_DIR=$TMP/${bam}.TMP
	printf "\n"

	echo "## BaseRecalibrator ..."
	time $GATK/gatk BaseRecalibrator \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $TMP/${bam}.TMP/${bam/.bam/.dedup.bam} \
		--known-sites $REF/mus_musculus.vcf \
		-O ../output/${bam/.bam/.dedup.bqsr.table} 
	printf "\n"

	echo "## ApplyBQSR ..."
	time $GATK/gatk ApplyBQSR \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $TMP/${bam}.TMP/${bam/.bam/.dedup.bam} \
		--bqsr-recal-file $RES/${bam/.bam/.dedup.bqsr.table} \
		-O ../output/${bam/.bam/.dedup.bqsr.bam}
	printf "\n"

	if [ -e $RES/${bam/.bam/.dedup.bqsr.bam} ]
	then
		echo "## Removing temporary files (keeping only ${bam/.bam/.dedup.bqsr.bam})"
		rm -r $TMP/${bam}.TMP
		printf "\n"

		echo "## Done: $bam converted into ${bam/.bam/.dedup.bqsr.bam}"
		printf "\n\n"
	else
		echo "## something went wrong while bqsr"
	fi
done

