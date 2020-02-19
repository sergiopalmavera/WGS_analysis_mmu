#!/bin/bash

# Set number of pairs 
FROM=1
TO=5

# Print version of programs
echo "# bwa version" # for bwa there is apparently no argument to query the version
bwa &> fl.tmp
grep Version fl.tmp
rm fl.tmp

echo "# Samtools version"
samtools --version

# Set number of threads
THR=20

# Set path to reference genome (relative path)
REF=../../../reference_genome_ensembl

# Set path to fastq "fastp-corrected" files
FASTQ=../../02_quality_trimming_adapter_removal/output

# Define name of R1 files
R1s=$(ls -1 $FASTQ | grep "fastq.gz" | grep "R1" | sed -n $FROM,${TO}p )

# Print a message with the list of pairs being processed
for R1 in $R1s; do echo "# Processing pairs ${R1/R1/[R1,R2]}"; done

# Loop over each R1 file
for R1 in $R1s
do
	# Print a message about the read pair being processed
	echo "## Processing pair ${R1/R1/[R1,R2]}"

	# Define R2 fastq file
	R2=$(echo ${R1/R1/R2})

	# Print a message with both fastq files in pair
	echo "### R1: $R1"
	echo "### R2: $R2"
	
	# Define name of BAM file	
	FLNM=$(echo $R1 | cut -d'_' -f1,2,3 | cut -d'.' -f1)
	
	# Run bwa-mem on fastq pair and pipe output to samtools to convert sam to bam
	echo "### Starting bwa mem"
	time bwa mem -t $THR -M $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz $FASTQ/$R1 $FASTQ/$R2 | samtools view -@ $THR -bS - > ../output/$FLNM.bam 
	printf "\n"

	echo "### Process completed, this is no confirmation that all went well, you need to inspect the BAM files!"
done
