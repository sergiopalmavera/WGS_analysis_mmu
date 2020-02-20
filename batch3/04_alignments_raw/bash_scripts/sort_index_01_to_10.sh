#!/bin/bash

# Define samples in batch
FROM=1
TO=10

# Absolute paths (modify accordingly)
PICARD=~/FBN_HOME/Tools/picard_2.18.11
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Relative paths
REF=../../../reference_genome_ensembl
RES=../output
TMP=../TMP

# Define samples
bam_files=$(ls -1 ../output/*.bam | sed -n $FROM,${TO}p)

echo "# Processing files"
for b in $bam_files; do echo $b; done
printf "\n"

time for b in $bam_files
do
	bam=$(basename $b) # define bam name
	echo "## BAM file $bam"
	
	mkdir $TMP/${bam}.TMP # make a temp dir for file
	
	echo "## Sorting ..."
	time java -jar $PICARD/picard.jar SortSam \
		I=$b \
		O=${b/.bam/.sorted.bam} \
		SORT_ORDER=coordinate \
		TMP_DIR=$TMP/${bam}.TMP
	printf "\n"

	echo "## Indexing ..."
	samtools index ${b/.bam/.sorted.bam}
	echo "\n"

	echo "## Removing temporary files and raw bam..."
	if [ -e ${b/.bam/.sorted.bam} ]
	then
		rm $b # remove raw bam
		rm -r $TMP/${bam}.TMP # remove temporary dir
		echo "## Done: $b converted into ${b/.bam/.sorted.bam} plus index ..."
		printf "\n\n"
	else
		echo "## $b was not sorted/idexed"
		printf "\n\n"
	fi

done
