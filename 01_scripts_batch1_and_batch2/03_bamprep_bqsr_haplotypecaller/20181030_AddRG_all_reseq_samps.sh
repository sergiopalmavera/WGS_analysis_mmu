#!/bin/bash

# Define paths to software
PICARD=~/FBN_HOME/Tools/picard_2.18.11
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

# Define path to reference genome
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

# Define path to fastqs, bams and a tmp dir multiple purpose
RAWFASTQ1=/projekte/I2-SOS-FERT/Original
RAWFASTQ2=/projekte/I2-SOS-FERT/Original2
OLDBAMS=/projekte/I2-SOS-FERT/03_alignments_raw/results
NEWBAMS=/projekte/I2-SOS-FERT/03_alignments_raw2/results
TMP=/projekte/I2-SOS-FERT/tmp

# Loop over each resequenced sample (reseq's due to low cvg)
cd $NEWBAMS
for file in *.bam
do
	# Define samples to process
	SAMPS=$(echo ${file%%-L*})

	# Collect BAMs for each sample
	for samp in $SAMPS
	do
		echo "#== Sample $samp ==#"
		echo "## Adding RG to BAMs for $samp (BAMs are temporary)"
		# Define bam files
		BAMorig=$(ls -1 $OLDBAMS | grep $samp)	
		BAMnew=$(ls -1 $NEWBAMS | grep $samp)
		
		# Make sample specific tmp dir
		mkdir $TMP/$samp
		mkdir $TMP/$samp/tmp

		# Add read groups to BAMs
		## original bams
		for i in $BAMorig
		do
			echo "## Adding RG to $i (BAM is temporary)"
		       	# Define FLOWCELL and LANE id
			FLOWCELL=$(zcat $RAWFASTQ1/${i/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d':' -f3)
			LANE=$(zcat $RAWFASTQ1/${i/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d ':' -f4)
			# Add RG to bam file
			time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
				I=$OLDBAMS/$i \
				O=$TMP/$samp/${i/bam/RG_BAMorig.bam} \
				RGID=$FLOWCELL.$LANE \
				RGLB=$(echo lib_$samp) \
				RGPL=$(echo ILLUMINA) \
				RGSM=$samp \
				RGPU=$FLOWCELL.$LANE.$samp
		done
		## new bams
		for i in $BAMnew
		do 
			echo "## Adding RG to $i (BAM is temporary)"
			# Define FLOWCELL and LANE id
			FLOWCELL=$(zcat $RAWFASTQ2/${i/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d':' -f3)
			LANE=$(zcat $RAWFASTQ2/${i/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d ':' -f4)
			# Add RG to bam file
			time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
				I=$NEWBAMS/$i \
				O=$TMP/$samp/${i/bam/RG_BAMnew.bam} \
				RGID=$FLOWCELL.$LANE \
				RGLB=$(echo lib_$samp) \
				RGPL=$(echo ILLUMINA) \
				RGSM=$samp \
				RGPU=$FLOWCELL.$LANE.$samp
		done
	done
done
