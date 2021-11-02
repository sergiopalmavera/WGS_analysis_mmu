#!/bin/bash

# Define date and print it
echo "Date of start"
date
DATE=`date +%Y%m%d_%H%M%S`

# Define path raw fastqs
RAWFASTQ=/projekte/I2-SOS-FERT/Original

# Define sample to process using original fastq directory as reference (consistent with former scripts)
FROM=11
TO=15
SAMPLES=$(ls -1 $RAWFASTQ | sed '/md5/d' | grep .fastq.gz | cut -d'_' -f1 | sort | uniq | sed -n $FROM,${TO}p)

# Define path to raw alignments
RAWBAM=/projekte/I2-SOS-FERT/03_alignments_raw/results

# Define path to output
RES=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/results

# Define path to tmp
TMP=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR/tmp

# Define path to PICARDTOOLS
PICARD=~/FBN_HOME/Tools/picard_2.18.11

# Define path to GATK
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Define path to reference genome and snps
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

# Define technology
PL=ILLUMINA

time for SAMP in $SAMPLES
do
	# Print sample name
	echo ===============================
	echo '## Sample:' $SAMP
	echo ==============================

	# Make a temporary directory for this sample, incl a tmp dir
	mkdir $TMP/$SAMP
	mkdir $TMP/$SAMP/tmp

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
			# add read groups
			echo "Adding read groups"
			time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
				I=$RAWBAM/$file \
				O=$TMP/$SAMP/${file%%.bam}.RG.bam \
				RGID=$FLOWCELL.$LANE \
				RGLB=$LB \
				RGPL=$PL \
				RGSM=$SAMP \
				RGPU=$FLOWCELL.$LANE.$SAMP
		done
		echo ==========================================
		echo "## Finished adding RG to BAMs for $SAMP"
		echo ==========================================

		# Collect bam names for merging
		SAMPBAMSPROC=$(echo $SAMPBAMS | sed 's/.bam/.RG.bam/g; s/'$SAMP'/I='$SAMP'/g')
		echo "Merging BAMs for $SAMP"

		# Change to processed-bams folder
		cd $TMP/$SAMP

		# Merge  
		echo "Merging"
		time java -jar $PICARD/picard.jar MergeSamFiles \
			$(echo $SAMPBAMSPROC) \
			O=$TMP/$SAMP/$SAMP.merged.bam \
			TMP_DIR=$TMP/$SAMP/tmp


		# remove temporary per lane bams 
		rm $TMP/$SAMP/*.RG.bam

		# Sort again
		echo "sorting bam"
		time java -jar $PICARD/picard.jar SortSam \
			I=$TMP/$SAMP/$SAMP.merged.bam \
			O=$TMP/$SAMP/$SAMP.merged.sorted.bam \
			SORT_ORDER=coordinate \
			TMP_DIR=$TMP/$SAMP/tmp

		# remove temporary file
		rm $TMP/$SAMP/$SAMP.merged.bam

		# Deduplicate merged bam 
		echo "Deduplicate"
		time java -jar $PICARD/picard.jar MarkDuplicates \
			 I=$TMP/$SAMP/$SAMP.merged.sorted.bam \
			 O=$TMP/$SAMP/$SAMP.merged.sorted.dedup.bam \
			 M=$RES/merged/$SAMP.merged.sorted.dedup.metrics.txt \
			 TMP_DIR=$TMP/$SAMP/tmp

		# Remove temporary merged bam
		rm $TMP/$SAMP/$SAMP.merged.sorted.bam

		# Base recalibration (BQSR)
		echo "Base recalibration: BaseRecalibrator"
		time $GATK/gatk BaseRecalibrator \
			-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
			-I $TMP/$SAMP/$SAMP.merged.sorted.dedup.bam \
			--known-sites $REF/mus_musculus.vcf \
			-O $TMP/$SAMP/$SAMP.merged.sorted.dedup.bqsr.table 

		# Move output to permamnent directory
		mv $TMP/$SAMP/$SAMP.merged.sorted.dedup.bqsr.table $RES/merged

		echo "Base reclaibration: ApplyBQSR"
		time $GATK/gatk ApplyBQSR \
			-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
			-I $TMP/$SAMP/$SAMP.merged.sorted.dedup.bam \
			--bqsr-recal-file $RES/merged/$SAMP.merged.sorted.dedup.bqsr.table \
			-O $TMP/$SAMP/$SAMP.merged.sorted.dedup.bqsr.bam

		# Move recalibrated bam into temporary directory
		mv $TMP/$SAMP/$SAMP.merged.sorted.dedup.bqsr.bam $RES/merged

		# remove temporary file
		rm $TMP/$SAMP/$SAMP.merged.sorted.dedup.bam
	else
		echo "The number of bam files for $SAMP is not 9"
	fi
	# Remove sample temporary file
	rm -r $TMP/$SAMP
done
