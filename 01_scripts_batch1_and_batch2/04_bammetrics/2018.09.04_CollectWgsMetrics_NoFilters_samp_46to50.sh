#!/bin/bash

# Define chunk
FROM=46
TO=50

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

# Set path to directory with BQSR-BAMs
BAM_BQSR=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/results/merged

# path to final results
OUT=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/CollectWgsMetrics_NoFilters/results

# Define path to tmp 
TMP=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/CollectWgsMetrics_NoFilters/tmp

# Define path to PICARDTOOLS
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

for file in $SAMPLES
do
	echo "## $file"
	printf "%s\n\n"
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$BAM_BQSR/$file.merged.sorted.dedup.bqsr.bam \
		O=$TMP/${file%%.bam}.CollectWgsMetrics_NoFilters.txt \
		R=$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		Q=0 \
	       	MQ=0	
	printf "%s\n\n\n\n"
	mv $TMP/${file%%.bam}.CollectWgsMetrics.txt $OUT
done


