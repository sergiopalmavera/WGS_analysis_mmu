#!/bin/bash

# Speciall treatment for this sample
# This belong to one of the resequenced samples
# One of the reseq BAM files has the same name as one of the original BAMs

# Define chunk
FROM=16
TO=16

# Define paths
BAMoriginals=/projekte/I2-SOS-FERT/03_alignments_raw/results
BAMreseq=/projekte/I2-SOS-FERT/03_alignments_raw2/results
PICARD=~/FBN_HOME/Tools/picard_2.18.11
TMP=/projekte/I2-SOS-FERT/tmp
OUT=/projekte/I2-SOS-FERT/12_alignments_raw_RG_merged_sorted_dedup

# Get sample names (60 samples)
RAWFASTQ=/projekte/I2-SOS-FERT/Original
cd $RAWFASTQ
SAMPLES=$(ls -1 | grep -v md5 | grep fastq.gz | for i in $(cat); do echo ${i:0:6}; done | sort | uniq)

# Define chunk
CHUNK=$(for i in $SAMPLES; do echo $i; done | sed -n $FROM,${TO}p)

echo "Looping through samples $FROM to $TO"
echo $CHUNK

# Loop over each sample
for SAMP in $CHUNK
do
	echo "# Processing bam files for sample $SAMP"
	
	# Collect bams for sample
	set1=$(ls -1 $BAMoriginals/*.bam | grep $SAMP)
	set2=$(ls -1 $BAMreseq/*.bam | grep $SAMP)
	SAMPBAMS=$(echo $set1 $set2)
	
	echo "## Number of raw bams for sample: $(echo $SAMPBAMS | wc -w)"

	# Add read groups to each bam file
	echo "## Adding read groups"
	for file in $SAMPBAMS
	do
		# Define output name avoiding duplicated names
		end=$(echo ${file%/*} | awk -F'/' '{print $4}' | awk '{print substr($0,length,1)}')
		nm=${file##*/}
		out_nm=$nm$end

		# Define flowcell and cell lane id
		FLOWCELL=$(zcat $RAWFASTQ/${nm/_001.bam/_R1_001.fastq.gz} | head -1 | cut -d ':' -f3)
		LANE=$(zcat $RAWFASTQ/${nm/_001.bam/_R1_001.fastq.gz}| head -1 | cut -d ':' -f4)

		# add read groups
		java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
			I=$file \
			O=$TMP/${out_nm/.bam/.RG.bam} \
			RGID=$FLOWCELL.$LANE \
			RGLB=lib_$SAMP \
			RGPL=ILLUMINA \
			RGSM=$SAMP \
			RGPU=$FLOWCELL.$LANE.$SAMP
	done
	
	cd $TMP

	echo "## Number of bams with RG:" 
	NRG=$(ls -1 | grep $SAMP | grep "RG")
	echo $NRG

	## Collect BAM files for sample
	SAMPBAMSPROC=$(ls -1 | grep $SAMP | sed 's/'$SAMP'/I='$SAMP'/g')

	echo "## Number of BAMs to merge"
	echo $SAMPBAMSPROC | wc -w

	## Make a tmp dir for this sample (remove at the end)
	mkdir $SAMP.tmp

	echo "## Merging and sorting BAMs with RG added"
	java -jar $PICARD/picard.jar MergeSamFiles \
		$(echo $SAMPBAMSPROC) \
		O=$SAMP.RG.merged.sorted.bam \
		SORT_ORDER=coordinate \
		TMP_DIR=$SAMP.tmp

	echo "## Marking duplicates"
	java -jar $PICARD/picard.jar MarkDuplicates \
		I=$SAMP.RG.merged.sorted.bam \
		O=$OUT/$SAMP.RG.merged.sorted.dedup.bam \
		M=$OUT/$SAMP.RG.merged.sorted.dedup.metrics \
		TMP_DIR=$SAMP.tmp

	## Remove temporary files
	if [ -e $SAMP.RG.merged.sorted.dedup.bam ]
	then
		rm -r $SAMP*
		echo "## Done processing sample $SAMP"
	fi
done
