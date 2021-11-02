#!/bin/bash
# Define date and print it
echo "Date of start"
date
DATE=`date +%Y%m%d_%H%M%S`

# Define path to raw alignments
RAWBAM=/projekte/I2-SOS-FERT/03_alignments_raw/results

# Define path raw fastqs
RAWFASTQ=/projekte/I2-SOS-FERT/Original

# Define path to output
RES=/projekte/I2-SOS-FERT/04_alignments_merged_dedup/results

# Define path to tmp
TMP=/projekte/I2-SOS-FERT/04_alignments_merged_dedup/tmp

# Define path to PICARDTOOLS
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Define path to GATK
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

# Define sample to process using original fastq directory as reference (consistent with former scripts)
FROM=1
TO=1
SAMPLES=$(ls -1 $RAWFASTQ | sed '/md5/d' | grep .fastq.gz | cut -d'_' -f1 | sort | uniq | sed -n $FROM,${TO}p)

# Define technology
PL=ILLUMINA

for SAMP in $SAMPLES
do
	# Print sample name
	echo ===============================
	echo '## Sample:' $SAMP
	echo ==============================

	# Define library name
	LB=lib_$SAMP

	# Search and count bam files for sample
	N=$(ls -1 $RAWBAM | grep -c $SAMP)

	# Capture bam names
	SAMPBAMS=$(ls -1 $RAWBAM | grep $SAMP)

	# Check if enough bams
	if [ $N = 9 ]
	then
		echo ========================================
		echo "## 9 BAM files for $SAMP found"
		echo ========================================
		for file in $SAMPBAMS
			do
			echo =======================================
			echo '## Sample ' $SAMP '; BAM ' $file
			echo =======================================
			# Define flow cell lane
			FLOWCELL=$(zcat $RAWFASTQ/${file/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d ':' -f3) 
			LANE=$(zcat $RAWFASTQ/${file/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d ':' -f4) 
			# Sort sam
			echo "sorting bam"
			time java -jar $PICARD/picard.jar SortSam \
				I=$RAWBAM/$file \
				O=$TMP/${file%%.bam}.sorted.bam \
				SORT_ORDER=coordinate
			# add read groups
			echo "Adding read groups"
			time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
				I=$TMP/${file%%.bam}.sorted.bam \
				O=$TMP/${file%%.bam}.sorted.RG.bam \
				RGID=$FLOWCELL.$LANE \
				RGLB=$LB \
				RGPL=$PL \
				RGSM=$SAMP \
				RGPU=$FLOWCELL.$LANE.$SAMP
			# Remove temporary bam
			rm $TMP/${file%%.bam}.sorted.bam
			# deduplicate 
			echo "Deduplicate by RG as QC"
			time java -jar $PICARD/picard.jar MarkDuplicates \
				I=$TMP/${file%%.bam}.sorted.RG.bam \
				O=$RES/single/${file%%.bam}.RG.sorted.dedup.bam \
				M=$RES/single/${file%%.bam}.RG.sorted.dedup.metrics.txt
			# remove temporary bam
			rm $TMP/${file%%.bam}.sorted.RG.bam
		done
		echo ==========================================
		echo "## Finished processing BAMs for $SAMP"
		echo ==========================================
	else
		echo "The number of bam files for $SAMP is not 9"
	fi

done

# Collect bam names for merging
SAMPBAMSPROC=$(echo $SAMPBAMS | sed 's/.bam/.RG.sorted.dedup.bam/g; s/'$SAMP'/I='$SAMP'/g')
echo "Merging BAMs for $SAMP"

# Change to processed-bams folder
cd $RES/single

# Merge  
echo "Merging"
time java -jar $PICARD/picard.jar MergeSamFiles \
	$(echo $SAMPBAMSPROC) \
	O=$TMP/$SAMP.merged.bam

# Sort again
echo "sorting bam"
time java -jar $PICARD/picard.jar SortSam \
	I=$TMP/$SAMP.merged.bam \
	O=$TMP/$SAMP.merged.sorted.bam \
	SORT_ORDER=coordinate

# remove tmp file
rm $TMP/$SAMP.merged.bam

# Deduplicate merged bam 
echo "Merging/deduplicate"
time java -jar $PICARD/picard.jar MarkDuplicates \
	 I=$TMP/$SAMP.merged.sorted.bam \
	 O=$TMP/$SAMP.merged.sorted.dedup.bam \
	 M=$TMP/$SAMP.merged.sorted.dedup.metrics.txt

# Remove temporary merged bam
rm $TMP/$SAMP.merged.sorted.bam

echo "Base reclaibration: ApplyBQSR"
time $GATK/gatk ApplyBQSR \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-I $TMP/$SAMP.merged.sorted.dedup.bam \
	--bqsr-recal-file $RES/merged/$SAMP.merged.sorted.dedup.recal.table \
        -O $RES/merged/$SAMP.merged.sorted.dedup.recal.bam

# Remove temporary file
rm $TMP/$SAMP.merged.sorted.dedup.bam #in case I need to rerun analysis


