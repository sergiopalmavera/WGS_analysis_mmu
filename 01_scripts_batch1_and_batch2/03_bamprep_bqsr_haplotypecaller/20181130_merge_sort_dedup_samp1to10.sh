#!/bin/bash

# Define chunk
FROM=1
TO=10

# Define paths
PICARD=~/FBN_HOME/Tools/picard_2.18.11
TMP=/projekte/I2-SOS-FERT/tmp
OUT=/projekte/I2-SOS-FERT/12_alignments_raw_RG_merged_sorted_dedup

# Get sample names (60 samples)
RAWFASTQ=/projekte/I2-SOS-FERT/Original
cd $RAWFASTQ
SAMPLES=$(ls -1 | grep -v md5 | grep fastq.gz | for i in $(cat); do echo ${i:0:6}; done | sort | uniq)

# Define chunk
CHUNK=$(for i in $SAMPLES; do echo $i; done | sed -n $FROM,${TO}p)

echo "Looping through samples $FROM to $TO"
echo $CHUNK

# Loop over each sample
for SAMP in $CHUNK
do
	echo "# Processing bam files for sample $SAMP"

	cd $TMP
		
	## Collect BAM files for sample	
	SAMPBAMSPROC=$(ls -1 | grep $SAMP | sed 's/'$SAMP'/I='$SAMP'/g')
	
	echo "## Number of BAMs to merge"
	echo $SAMPBAMSPROC | wc -w
	
	## Make a tmp dir for this sample (remove at the end)
	mkdir $SAMP.tmp

	echo "## Merging and BAMs with RG added"
	java -jar $PICARD/picard.jar MergeSamFiles \
		$(echo $SAMPBAMSPROC) \
		O=$SAMP.RG.merged.sorted.bam \
		SORT_ORDER=coordinate \
		TMP_DIR=$SAMP.tmp

	echo "## Marking duplicates"
	java -jar $PICARD/picard.jar MarkDuplicates \
	       I=$SAMP.RG.merged.sorted.bam \
	       O=$OUT/$SAMP.RG.merged.sorted.dedup.bam \
	       M=$OUT/$SAMP.RG.merged.sorted.dedup.metrics \
	       TMP_DIR=$SAMP.tmp
	
	## Remove temporary files
	if [ -e $SAMP.RG.merged.sorted.dedup.bam ]
	then
		rm -r $SAMP*
		echo "## Done processing sample $SAMP"
	fi
done
