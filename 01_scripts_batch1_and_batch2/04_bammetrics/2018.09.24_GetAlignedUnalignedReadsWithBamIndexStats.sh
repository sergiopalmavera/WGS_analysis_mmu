PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
BAMs=/projekte/I2-SOS-FERT/03_alignments_raw/merged
OUT=/projekte/I2-SOS-FERT/03_alignments_raw/merged_BamIndexStats_flagstat
for bam in $(ls -1 $BAMs/*.bam)
do
	NM=${bam##*/}
	echo "##" $NM
	echo "## Starting samtools index"
	time samtools index $bam
	echo "## starting samtools flagstat"
	time samtools flagstat $bam > $OUT/${NM%.bam}.flagstat.tab
	echo "## starting picard tools BamIndexStats"
	time java -jar $PICARD/picard.jar BamIndexStats INPUT=$bam > $OUT/${NM%.bam}.BAMIndexStats.tab
done
