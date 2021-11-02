REF=/projekte/I2-SOS-FERT/reference_genome_ensembl
PICARD=~/FBN_HOME/Tools/picard_2.18.11
RES_FILTERS=/projekte/I2-SOS-FERT/03_alignments_raw2/CollectWgsMetrics_with_filters
RES_NO_FILTERS=/projekte/I2-SOS-FERT/03_alignments_raw2/CollectWgsMetrics_no_filters

cd /projekte/I2-SOS-FERT/03_alignments_raw2/merged
for bam in *.bam
do
	echo "# Calculating CollectWgsMetrics on $bam"
	echo "## CollectWgsMetrics with default (Q=20 & MQ=20) filters for $bam"
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$bam \
		O=$RES_FILTERS/${bam%%.bam}.CollectWgsMetrics_WithFilters.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa

	echo "## CollectWgsMetrics without filters for $bam"
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$bam \
		O=$RES_NO_FILTERS/${bam%%.bam}.CollectWgsMetrics_NoFilters.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		Q=0 \
		MQ=0
done
