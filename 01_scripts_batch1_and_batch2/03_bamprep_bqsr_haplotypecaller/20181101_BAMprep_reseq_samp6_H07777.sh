# Define paths
TMP=/projekte/I2-SOS-FERT/tmp
RES=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR2/results
MERGED=/projekte/I2-SOS-FERT/03_alignments_raw2/merged
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Define sample:
samp=H07777

echo "# Working on $samp"
# Prepend I= to bam files for merging
BAMs=$(for i in $(ls -1 $TMP/$samp/*.bam); do echo "I=$i"; done)

# Make a samp tmp for picard
mkdir $TMP/tmp/$samp

# Merge and remove inputs if succesfull
echo "## Merging BAMs for $samp"
time java -jar $PICARD/picard.jar MergeSamFiles \
	$(echo $BAMs) \
	O=$TMP/$samp.merged.bam \
	TMP_DIR=$TMP/tmp/$samp

if [ -e $TMP/$samp.merged.bam ]; then rm -r $TMP/$samp; fi

# Sort and then move input if succesful
echo "## Sorting merged BAM for $samp"
time java -jar $PICARD/picard.jar SortSam \
	I=$TMP/$samp.merged.bam \
	O=$TMP/$samp.merged.sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=$TMP/tmp/$samp

if [ -e $TMP/$samp.merged.sorted.bam ]; then mv $TMP/$samp.merged.bam $MERGED; fi

# Deduplicate and remove input if succesful
echo "## Deduplicating BAM for $samp"
time java -jar $PICARD/picard.jar MarkDuplicates \
	I=$TMP/$samp.merged.sorted.bam \
	O=$TMP/$samp.merged.sorted.dedup.bam \
	M=$RES/$samp.merged.sorted.dedup.metrics.txt \
	TMP_DIR=$TMP/tmp/$samp

if [ -e $TMP/$samp.merged.sorted.dedup.bam ]
then 
	echo "## BAM file for $samp ready for BQSR"
	rm $TMP/$samp.merged.sorted.bam #also remove input file
	rm -r $TMP/tmp/$samp #remove picard's sample tmp
fi
