#!/bin/bash

# Define paths 
TMP=/projekte/I2-SOS-FERT/tmp
BQSRtable=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR2/BQSRtable
BQSRbam=/projekte/I2-SOS-FERT/04_alignments_merged_dedup_BQSR2/BQSRbam
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

cd $TMP
time for SAMP in *.bam
do
	echo "# Running BQSR on $SAMP"
	echo "## Part1: BaseRecalibrator"
	time $GATK/gatk BaseRecalibrator \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $SAMP \
		--known-sites $REF/mus_musculus.vcf \
		-O $BQSRtable/${SAMP%%.bam}.bqsr.table 

	echo "## Part2: ApplyBQSR on $SAMP"
	time $GATK/gatk ApplyBQSR \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-I $SAMP \
		--bqsr-recal-file $BQSRtable/${SAMP%%.bam}.bqsr.table \
		-O $BQSRbam/${SAMP%%.bam}.bqsr.bam

	if [ -e $BQSRbam/${SAMP%%.bam}.bqsr.bam ]
	then 
		echo "## BQSR for $SAMP done"
		echo "## $SAMP has been removed to reduce storage"
		rm $SAMP
	fi 
done
