#!/bin/bash

# Absolute paths (modify accordingly)
PICARD=~/FBN_HOME/Tools/picard_2.18.11
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
FASTQ=/projekte/I2-SOS-FERT/Original3

# Relative paths
bam_files_path=../../04_alignments_raw/output
REF=../../../reference_genome_ensembl
TMP=../TMP

bam=I34772-L1_S19_L004.sorted.bam
b=$bam_files_path/I34772-L1_S19_L004.sorted.bam
	
echo "## BAM file $bam"
mkdir $TMP/${bam}.TMP

echo "## Adding Read Group ..."
fastq=${bam/.sorted.bam/_R1_001.fastq.gz} # recreate fastq file name
FLOWCELL=$(zcat $FASTQ/$fastq | head -1 | cut -d ':' -f3)
LANE=$(zcat $FASTQ/$fastq | head -1 | cut -d ':' -f4)
samp=$( echo $bam | cut -d'-' -f1 )
time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
	I=$b \
	O=$TMP/${bam}.TMP/${bam/.bam/.RG.bam} \
	RGID=$FLOWCELL.$LANE \
	RGLB=lib_$samp \
	RGPL=ILLUMINA \
	RGSM=$samp \
	RGPU=$FLOWCELL.$LANE.$samp

echo "## Marking duplicates ..."
time java -jar $PICARD/picard.jar MarkDuplicates \
	I=$TMP/${bam}.TMP/${bam/.bam/.RG.bam} \
	O=$TMP/${bam}.TMP/${bam/.bam/.RG.dedup.bam} \
	M=../output/${bam/.bam/.dedup.metrics.txt} \
	TMP_DIR=$TMP/${bam}.TMP
printf "\n"

echo "## BaseRecalibrator ..."
time $GATK/gatk BaseRecalibrator \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-I $TMP/${bam}.TMP/${bam/.bam/.RG.dedup.bam} \
	--known-sites $REF/mus_musculus.vcf \
	-O ../output/${bam/.bam/.RG.dedup.bqsr.table} 
printf "\n"

echo "## ApplyBQSR ..."
time $GATK/gatk ApplyBQSR \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-I $TMP/${bam}.TMP/${bam/.bam/.RG.dedup.bam} \
	--bqsr-recal-file ../output/${bam/.bam/.RG.dedup.bqsr.table} \
	-O ../output/${bam/.bam/.RG.dedup.bqsr.bam}
printf "\n"

if [ -e ../output/${bam/.bam/.RG.dedup.bqsr.bam} ]
then
	echo "## Done: $bam converted into ${bam/.bam/.RG.dedup.bqsr.bam}"
	printf "\n\n"
else
	echo "## something went wrong while bqsr"
	printf "\n\n"
fi

rm -r $TMP/${bam}.TMP # remove temp files

