#!/bin/bash

# Define chunk
FROM=1
TO=10

# Print hostname
echo "Hostname"
hostname

# Define date and print it
echo "Date of start"
date
DATE=`date +%Y%m%d_%H%M%S`

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
OUT_CollectWgsMetrics_WithFilters=/projekte/I2-SOS-FERT/03_alignments_raw/CollectWgsMetrix_WithFilters/results
OUT_SamtoolsDepth_WithFilters=/projekte/I2-SOS-FERT/03_alignments_raw/DepthAlongGenome_WithFilters/results
OUT_CollectWgsMetrics_NoFilters=/projekte/I2-SOS-FERT/03_alignments_raw/CollectWgsMetrix_NoFilters/results
OUT_SamtoolsDepth_NoFilters=/projekte/I2-SOS-FERT/03_alignments_raw/DepthAlongGenome_NoFilters/results

# Define path to merged.sorted.dedup BAMs
MERGED=/projekte/I2-SOS-FERT/03_alignments_raw/merged

for SAMP in $SAMPLES
do
	# Make a temporary directory for this sample, incl a tmp dir
	mkdir $TMP/$SAMP
	mkdir $TMP/$SAMP/tmp

	# Get BAMs for this sample
	SAMPBAMS=$(ls -1 $RAWBAM | grep $SAMP)
		
 	# Merge BAMs
	java -jar $PICARD/picard.jar MergeSamFiles \
		$(echo $SAMPBAMS) \
		O=$TMP/$SAMP/$SAMP.merged.bam \
		TMP_DIR=$TMP/$SAMP/tmp

	# Sort
	java -jar $PICARD/picard.jar SortSam \
		I=$TMP/$SAMP/$SAMP.merged.bam \
		O=$TMP/$SAMP/$SAMP.merged.sorted.bam \
		SORT_ORDER=coordinate \
		TMP_DIR=$TMP/$SAMP/tmp
	rm $TMP/$SAMP/$SAMP.merged.bam

	# Mark duplicates
	java -jar $PICARD/picard.jar MarkDuplicates \
		 I=$TMP/$SAMP/$SAMP.merged.sorted.bam \
		 O=$MERGED/$SAMP.merged.sorted.dedup.bam \
		 M=$MERGED/$SAMP.merged.sorted.dedup.metrics.txt \
		 TMP_DIR=$TMP/$SAMP/tmp
	rm $TMP/$SAMP/$SAMP.merged.sorted.bam

	# CollectWgsMetrics with filters
	java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$MERGED/$SAMP.merged.sorted.dedup.bam \
		O=$OUT_CollectWgsMetrics_WithFilters/$SAMP.merged.CollectWgsMetrics_WithFilters.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa

	# Samtools depth with filters
	samtools view -b -F 0x400 $MERGED/$SAMP.merged.sorted.dedup.bam > $TMP/$SAMP/$SAMP.tmp.bam 
	samtools depth -q 20 -Q 20 $TMP/$SAMP/$SAMP.tmp.bam > $OUT_SamtoolsDepth_WithFilters/$SAMP.merged.sorted.dedup_WithFilters.cvg #Qq same as picard's 

	# CollectWgsMetrics No filters
	java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$MERGED/$SAMP.merged.sorted.dedup.bam \
		O=$OUT_CollectWgsMetrics_NoFilters/$SAMP.merged.CollectWgsMetrics_NoFilters.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		Q=0 \
		MQ=0

	# Samtools depth No filters
	samtools depth -q 0 -Q 0 $TMP/$SAMP/$SAMP.tmp.bam > $OUT_SamtoolsDepth_NoFilters/$SAMP.merged.sorted.dedup_NoFilters.cvg #Qq same as picard's 
	rm $TMP/$SAMP/$SAMP.tmp.bam

	# remove temporary files
	rm -r $TMP/$SAMP
done
