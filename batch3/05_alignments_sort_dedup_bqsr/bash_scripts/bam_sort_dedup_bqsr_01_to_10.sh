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
RES=../output
TMP=../TMP

# Define samples
bam_files=$(ls -1 $bam_files_path/*.bam | sed -n $FROM,${TO}p)

echo "# Processing files"
for bam in $bam_files; do echo $bam; done
printf "\n"

time for b in $bam_files
do
	bam=$(basename $b)
	
	echo "## BAM file $bam"
	mkdir $TMP/${bam}.TMP
	
	echo "## Sorting ..."
	time java -jar $PICARD/picard.jar SortSam \
		I=$b \
		O=$TMP/${bam}.TMP/${bam/.bam/.sorted.bam} \
		SORT_ORDER=coordinate \
		TMP_DIR=$TMP/${bam}.TMP
	printf "\n"

	echo "## Marking duplicates ..."
	time java -jar $PICARD/picard.jar MarkDuplicates \
		I=$TMP/${bam}.TMP/${bam/.bam/.sorted.bam} \
		O=$TMP/${bam}.TMP/${bam/.bam/.sorted.dedup.bam} \
		M=$RES/${bam/.bam/.sorted.dedup.metrics.txt} \
		TMP_DIR=$TMP/${bam}.TMP
	printf "\n"

	echo "## BaseRecalibrator ..."
	time $GATK/gatk BaseRecalibrator \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $TMP/${bam}.TMP/${bam/.bam/.sorted.dedup.bam} \
		--known-sites $REF/mus_musculus.vcf \
		-O $RES/${bam/.bam/.sorted.dedup.bqsr.table} 
	printf "\n"

	echo "## ApplyBQSR ..."
	time $GATK/gatk ApplyBQSR \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $TMP/${bam}.TMP/${bam/.bam/.sorted.dedup.bam} \
		--bqsr-recal-file $RES/${bam/.bam/.sorted.dedup.bqsr.table} \
		-O $RES/${bam/.bam/.sorted.dedup.bqsr.bam}
	printf "\n"

	echo "## Removing temporary files (keeping only ${bam/.bam/.sorted.dedup.bqsr.bam})"
	rm -r $TMP/${bam}.TMP
	printf "\n"

	echo "## Done: $bam converted into ${bam/.bam/.sorted.dedup.bqsr.bam}"
	printf "\n\n"
done

