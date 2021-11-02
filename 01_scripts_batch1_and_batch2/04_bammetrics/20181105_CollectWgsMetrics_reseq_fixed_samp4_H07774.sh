# Define paths
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl
PICARD=~/FBN_HOME/Tools/picard_2.18.11
RES_FILTERS=/projekte/I2-SOS-FERT/03_alignments_raw2/CollectWgsMetrics_with_filters
RES_NO_FILTERS=/projekte/I2-SOS-FERT/03_alignments_raw2/CollectWgsMetrics_no_filters
TMP=/projekte/I2-SOS-FERT/tmp
MERGED_DEDUP=/projekte/I2-SOS-FERT/03_alignments_raw2/merged

# Define sample and corresponding bam

SAMP=H07774

cd /projekte/I2-SOS-FERT/03_alignments_raw2/merged
bam=$(ls -1 | grep $SAMP)

echo "# Processing merged bam file for sample $SAMP"
echo "## merged bam file $bam"

# Make a samp tmp for picard
mkdir $TMP/tmp/$SAMP

# Sort and then remove input if succesful
echo "## Sorting merged BAM for $samp"
time java -jar $PICARD/picard.jar SortSam \
	I=$bam \
	O=$TMP/${bam%%.bam}.sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=$TMP/tmp/$samp

if [ -e $TMP/${bam%%.bam}.sorted.bam ]; then rm $bam; fi

# Mark duplicates (this step was omited before and now is added, thus fixing the problem)
echo "## Deduplicating BAM for $samp"
time java -jar $PICARD/picard.jar MarkDuplicates \
	I=$TMP/${bam%%.bam}.sorted.bam \
	O=$TMP/${bam%%.bam}.sorted.dedup.bam \
	M=$TMP/$SAMP.metrics.txt \
	TMP_DIR=$TMP/tmp/$samp

if [ -e $TMP/${bam%%.bam}.sorted.dedup.bam ]; then rm $TMP/${bam%%.bam}.sorted.bam; fi
if [ -e $TMP/${bam%%.bam}.sorted.dedup.bam ]; then rm $TMP/$SAMP.metrics.txt; fi #metrics were produced earlier

# Apply CollectWgsMetrics 
echo "# Calculating CollectWgsMetrics on ${bam%%.bam}.sorted.dedup.bam"
echo "## CollectWgsMetrics with default (Q=20 & MQ=20) filters for $bam"
time java -jar $PICARD/picard.jar CollectWgsMetrics \
	I=$TMP/${bam%%.bam}.sorted.dedup.bam \
	O=$RES_FILTERS/${bam%%.bam}.sorted.dedup.CollectWgsMetrics_WithFilters.txt \
	R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa

echo "## CollectWgsMetrics without filters for $bam"
time java -jar $PICARD/picard.jar CollectWgsMetrics \
	I=$TMP/${bam%%.bam}.sorted.dedup.bam \
	O=$RES_NO_FILTERS/${bam%%.bam}.sorted.dedup.CollectWgsMetrics_NoFilters.txt \
	R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	Q=0 \
	MQ=0

if [[ -e $RES_FILTERS/${bam%%.bam}.sorted.dedup.CollectWgsMetrics_WithFilters.txt && -e $RES_FILTERS/${bam%%.bam}.sorted.dedup.CollectWgsMetrics_NoFilters.txt ]]; then mv $TMP/${bam%%.bam}.sorted.dedup.bam $MERGED_DEDUP; fi

rm -r $TMP/tmp/$SAMP 
