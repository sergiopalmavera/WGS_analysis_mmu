PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
BAMs=/projekte/I2-SOS-FERT/03_alignments_raw/merged
OUT=/projekte/I2-SOS-FERT/03_alignments_raw/merged_CollectInsertSizeMetrics

for bam in $(ls -1 $BAMs/*.bam)
do
	echo "##" $bam
	NM=${bam##*/}
	java -jar $PICARD/picard.jar CollectInsertSizeMetrics \
		I=$bam \
		O=$OUT/${NM%.bam}.CollectInsertSizeMetrics.tab \
		H=$OUT/${NM%.bam}.insert_size_hist.pdf
done

