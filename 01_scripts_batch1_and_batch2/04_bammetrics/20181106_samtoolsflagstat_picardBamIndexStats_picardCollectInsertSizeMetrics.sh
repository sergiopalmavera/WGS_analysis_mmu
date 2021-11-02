PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
BAMs=/projekte/I2-SOS-FERT/03_alignments_raw2/merged
OUT=/projekte/I2-SOS-FERT/03_alignments_raw2
for bam in $(ls -1 $BAMs/*.bam)
do
	# Extract name of bam file without path
	NM=${bam##*/}
	echo "# Processing BAM $NM"

	echo "## Running samtools index"
	time samtools index $bam

	echo "## Running samtools flagstat"
	time samtools flagstat $bam > $OUT/merged_BamIndexStats_flagstat/${NM%.bam}.flagstat.tab

	echo "## starting picard tools BamIndexStats"
	time java -jar $PICARD/picard.jar BamIndexStats INPUT=$bam > $OUT/merged_BamIndexStats_flagstat/${NM%.bam}.BAMIndexStats.tab

	echo "## Running Picard CollectInsertSizeMetrics"
	time java -jar $PICARD/picard.jar CollectInsertSizeMetrics \
		I=$bam \
		O=$OUT/merged_CollectInsertSizeMetrics/${NM%.bam}.CollectInsertSizeMetrics.tab \
		H=$OUT/merged_CollectInsertSizeMetrics/${NM%.bam}.insert_size_hist.pdf
done
