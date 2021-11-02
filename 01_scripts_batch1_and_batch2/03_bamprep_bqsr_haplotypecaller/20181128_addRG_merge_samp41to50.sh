#!/bin/bash

# Define chunk
FROM=41
TO=50

# Define paths
BAMoriginals=/projekte/I2-SOS-FERT/03_alignments_raw/results
BAMreseq=/projekte/I2-SOS-FERT/03_alignments_raw2/results
PICARD=~/FBN_HOME/Tools/picard_2.18.11
TMP=/projekte/I2-SOS-FERT/tmp

# Get sample names (60 samples)
RAWFASTQ=/projekte/I2-SOS-FERT/Original
cd $RAWFASTQ
SAMPLES=$(ls -1 | grep -v md5 | grep fastq.gz | for i in $(cat); do echo ${i:0:6}; done | sort | uniq)

# Define chunk
CHUNK=$(for i in $SAMPLES; do echo $i; done | sed -n $FROM,${TO}p)

echo "Looping through samples $FROM to $TO"
echo $CHUNK

# Loop over each sample
for samp in $CHUNK
do
	echo "# Processing bam files for sample $samp"
	
	# Collect bams for sample
	set1=$(ls -1 $BAMoriginals/*.bam | grep $samp)
	set2=$(ls -1 $BAMreseq/*.bam | grep $samp)
	SAMPBAMS=$(echo $set1 $set2)
	
	echo "## Number of bams for sample: $(echo $SAMPBAMS | wc -w)"

	# Add read groups to each bam file
	echo "## Adding read groups"
	for file in $SAMPBAMS
	do
		# Define flowcell and cell lane id
		nm=${file##*/}
		FLOWCELL=$(zcat $RAWFASTQ/${nm/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d ':' -f3)
		LANE=$(zcat $RAWFASTQ/${nm/_001.bam/_R1_001.fastq.gz}| head -1 | cut -d ':' -f4)

		# add read groups
		java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
			I=$file \
			O=$TMP/${nm/.bam/.RG.bam} \
			RGID=$FLOWCELL.$LANE \
			RGLB=lib_$samp \
			RGPL=ILLUMINA \
			RGSM=$samp \
			RGPU=$FLOWCELL.$LANE.$samp
	done

	echo "## Done processing sample $samp"
done
