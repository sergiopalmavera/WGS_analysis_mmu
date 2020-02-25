#!/bin/bash

FROM=11
TO=20

echo "# Working on:"
hostname
printf "\n"

# Retrieving metrics from bams sorted-AddRG-dedup-bqsr

# Absoulute path (modify accordingly)
picard=~/FBN_HOME/Tools/picard_2.18.11

bams=$(ls -1 ../../05_alignments_addRG_dedup_bqsr/output/*.bam | sed -n $FROM,${TO}p)
echo "# Processing samples:"
for b in $bams; do echo $b; done
printf "\n"

for b in $bams
do
	bam=$(basename $b)
	echo "## Processing $bam"
	printf "\n"

	echo "### samtools: flagstat ..."
	time samtools flagstat $b > ../output/${bam/.bam/.flagstat.tab}
	printf "\n"

	echo "### picard: BamIndexStats ..."
	time java -jar $picard/picard.jar BamIndexStats I=$b > ../output/${bam/.bam/.BAMIndexStats.tab}
	printf "\n"

	echo "### CollectWgsMetrics default ..."
	java -jar $picard/picard.jar CollectWgsMetrics \
		I=$b \
		O=../output/${bam/.bam/.CollectWgsMetrics.default.txt} \
		R=../../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa
	printf "\n"

	echo "### CollectWgsMetrics no filters..."
	java -jar $picard/picard.jar CollectWgsMetrics \
		I=$b \
		O=../output/${bam/.bam/.CollectWgsMetrics.Q0.M0.txt} \
		R=../../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		Q=0 \
		MQ=0
	printf "\n"
	
	echo "### picard: CollectInsertSizeMetrics ..."
	java -jar $picard/picard.jar CollectInsertSizeMetrics \
		I=$b \
		O=../output/${bam/.bam/.CollectInsertSizeMetrics.txt} \
		H=../output/${bam/.bam/.CollectInsertSizeMetrics.pdf}
	printf "\n"

done
