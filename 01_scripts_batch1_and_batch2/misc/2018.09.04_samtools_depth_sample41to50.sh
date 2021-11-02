#!/bin/bash

# Define chunk
FROM=41
TO=50

# Print hostname
echo "Hostname"
hostname

# Define date and print it
echo "Date of start"
date

# Path to BAM files
BAMDIR=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/results/merged

# Path to temporary directory
TMP=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/DepthAlongGenome/tmp

# Path to final directory
RES=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/DepthAlongGenome/results

# Define path raw fastqs
RAWFASTQ=/projekte/I2-SOS-FERT/Original

# Define sample to process using original fastq directory as reference (consistent with former scripts)
SAMPLES=$(ls -1 $RAWFASTQ | sed '/md5/d' | grep .fastq.gz | cut -d'_' -f1 | sort | uniq | sed -n $FROM,${TO}p)

time for samp in $SAMPLES
do
	echo "## SAMPLE $samp"	
	time samtools view -b -F 0x400 $BAMDIR/$samp.merged.sorted.dedup.bqsr.bam > $TMP/$samp.tmp.bam 
	time samtools depth -q 20 -Q 20 $TMP/$samp.tmp.bam > $TMP/$samp.merged.sorted.dedup.bqsr.cvg #Qq same as picard's CollectWgsMetrics
	rm $TMP/$samp.tmp.bam
	mv $TMP/$samp.merged.sorted.dedup.bqsr.cvg $RES
done
