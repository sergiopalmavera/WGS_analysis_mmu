
# Define paths
TMP=/projekte/I2-SOS-FERT/tmp
RES=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR2/results
MERGED=/projekte/I2-SOS-FERT/03_alignments_raw2/merged
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Get samples to process
samp=H07739

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


