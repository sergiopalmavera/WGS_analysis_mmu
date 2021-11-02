#!/bin/bash

# Define start and end of chunk
FROM=1
TO=10

## Defining paths
TMP=/projekte/I2-SOS-FERT/03_alignments_raw2/tmp
TRIMMD=/projekte/I2-SOS-FERT/02_trimmed2/results
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl
RES=/projekte/I2-SOS-FERT/03_alignments_raw2/results

## Set noumber of threads
THR=64

## Defining R1 files for chunk
ls -1 $TRIMMD | grep "OutputPaired" | grep "R1" | sed -n $FROM,${TO}p > ~/tmp/R1.reseq.bwa.$FROM.$TO

## Print out R1 one files
cat ~/tmp/R1.reseq.bwa.$FROM.$TO

while read file
do
	echo "#== Processing pair =="
	R1=$file
	R2=${file/R1/R2}
	echo $R1
	echo $R2

	FLNM=$(echo ${R1%%.Output*} | sed 's/R1_//g')
	
	echo "Looping over R1 and R2 pairs making first sure both files exist"
	if [ -e $TRIMMD/$R1 ] && [ -e $TRIMMD/$R2 ]
	then
		echo "## Both files exist"
		echo "## Runing BWA"
		time bwa mem -t $THR -M $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz $TRIMMD/$R1 $TRIMMD/$R2 > $TMP/$FLNM.sam
		echo "## Converting sam to bam with samtools"
		samtools --version
		time samtools view -@ $THR -bS $TMP/$FLNM.sam > $RES/$FLNM.bam
		echo "## Clearing tmp"
		rm $TMP/*
	else
		"Pair not found for $file"
	fi
done < ~/tmp/R1.reseq.bwa.$FROM.$TO
echo "#== End of chunk for pairs $FROM to $TO ==#"
