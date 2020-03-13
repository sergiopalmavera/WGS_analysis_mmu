#!/bin/bash

# Absolute paths (modify accordingly)
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Relative paths
REF=../../../reference_genome_ensembl
TMP=../TMP
OUT=../output

# Define samples
b=I34772-L1_S19_L004.bam

mkdir $TMP/${b}.TMP # make a temp dir for file

echo "## Sorting ..."
time java -jar $PICARD/picard.jar SortSam \
	I=$OUT/$b \
	O=$OUT/${b/.bam/.sorted.bam} \
	SORT_ORDER=coordinate \
	TMP_DIR=$TMP/${b}.TMP
printf "\n"

echo "## Indexing ..."
samtools index ../output/${b/.bam/.sorted.bam}
echo "\n"

rm $b # remove raw bam
rm -r $TMP/${b}.TMP # remove temporary dir
