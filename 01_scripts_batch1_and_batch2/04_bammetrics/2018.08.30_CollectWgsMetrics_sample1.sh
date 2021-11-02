#!/bin/bash

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
FROM=1
TO=1
SAMPLES=$(ls -1 $RAWFASTQ | sed '/md5/d' | grep .fastq.gz | cut -d'_' -f1 | sort | uniq | sed -n $FROM,${TO}p)

# Set path to directory with BQSR-BAMs
BAM_BQSR=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/results/merged

# path to final results
OUT=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/CollectWgsMetrics/results

# Define path to tmp 
TMP=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/CollectWgsMetrics/tmp

# Define path to PICARDTOOLS
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

for file in $SAMPLES
do
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$BAM_BQSR/$file.merged.sorted.dedup.bqsr.bam \
		O=$TMP/${file%%.bam}.CollectWgsMetrics.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa
	mv $TMP/${file%%.bam}.CollectWgsMetrics.txt $OUT
done


