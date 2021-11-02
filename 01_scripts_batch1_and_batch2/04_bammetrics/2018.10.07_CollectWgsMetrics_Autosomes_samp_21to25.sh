#!/bin/bash

# Define chunk
FROM=21
TO=25

# Define path raw fastqs
RAWFASTQ=/projekte/I2-SOS-FERT/Original

# Define sample to process using original fastq directory as reference (consistent with former scripts)
SAMPLES=$(ls -1 $RAWFASTQ | sed '/md5/d' | grep .fastq.gz | cut -d'_' -f1 | sort | uniq | sed -n $FROM,${TO}p)

# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

# Define path to PICARDTOOLS
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Define path to raw alignments
RAWBAM=/projekte/I2-SOS-FERT/03_alignments_raw/results

# Define path to TMP
TMP=/projekte/I2-SOS-FERT/03_alignments_raw/tmp

# Define output directories
OUT_CollectWgsMetrics_WithFilters=/projekte/I2-SOS-FERT/03_alignments_raw/CollectWgsMetrix_WithFilters_Autosomes/results
OUT_CollectWgsMetrics_NoFilters=/projekte/I2-SOS-FERT/03_alignments_raw/CollectWgsMetrix_NoFilters_Autosomes/results

# Define path to merged.sorted.dedup BAMs
MERGED=/projekte/I2-SOS-FERT/03_alignments_raw/merged

time for SAMP in $SAMPLES
do
	echo "# $SAMP"

	# Make a temp dir to store autosome-bam
	mkdir $TMP/$SAMP
	
	# Select only autosomes	
	echo "## Removing non-autosomes"
	time samtools idxstats $MERGED/$SAMP.merged.sorted.dedup.bam | cut -f 1 | grep -E -v [A-Za-z] | grep -E [0-9] | xargs samtools view -b $MERGED/$SAMP.merged.sorted.dedup.bam >$TMP/$SAMP/$SAMP.tmp.bam 	

	# CollectWgsMetrics with filters
	echo "## Collecting metrics (baseQ 20 and MAPQ 20"
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$TMP/$SAMP/$SAMP.tmp.bam \
		O=$OUT_CollectWgsMetrics_WithFilters/$SAMP.merged.sorted.dedup.CollectWgsMetrics_WithFilters_Autosomes.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa

	# CollectWgsMetrics No filters
	echo "## Collectiong metrics (baseQ 0 MAPQ 0)"
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$TMP/$SAMP/$SAMP.tmp.bam \
		O=$OUT_CollectWgsMetrics_NoFilters/$SAMP.merged.sorted.dedup.CollectWgsMetrics_NoFilters_Autosomes.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		Q=0 \
		MQ=0

	# Remove temporary BAM	
	rm $TMP/$SAMP/$SAMP.tmp.bam

	# remove temporary sample directory
	rm -r $TMP/$SAMP

	# Print finishing message
	echo "## Done calculating metrics for" $SAMP
done
