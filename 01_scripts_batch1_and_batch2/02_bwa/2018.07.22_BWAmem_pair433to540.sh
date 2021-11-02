#!/bin/bash

# Define start and end of chunk
FROM=433
TO=540

# Define a counter variable to print out the number of iterations
count=$FROM

echo "Start of chunk for pairs $FROM to $TO"
echo "Defining paths"
TMP=/projekte/I2-SOS-FERT/03_alignments_raw/tmp
TRIMMD=/projekte/I2-SOS-FERT/02_trimmed/results/TRIMMO
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

THR=64

echo "Defining R1 files for chunk"
ls -1 $TRIMMD | grep "OutputPaired" | grep "R1" | sed -n $FROM,${TO}p > ~/tmp/R1.bwa.$FROM.$TO

echo "Date and time starting of chunk pairs $FROM to $TO out of 540 pairs"
date

echo "Creating temporary directory for chunk"
mkdir -p $TMP/tmp.$FROM.$TO

while read file
do
	echo "#== Processing pair $count =="
	R1=$file
	R2=${file/R1/R2}
	echo $R1
	echo $R2

	FLNM=$(echo ${R1%%.Output*} | sed 's/R1_//g')
	
	echo "Looping over R1 and R2 pairs making first sure both files exist"
	if [ -e $TRIMMD/$R1 ] && [ -e $TRIMMD/$R2 ]
	then
		echo "Both files exist"
		echo "Run BWA"
		time bwa mem -t $THR -M $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz $TRIMMD/$R1 $TRIMMD/$R2 > $TMP/tmp.$FROM.$TO/$FLNM.sam
		echo "Converting sam to bam with samtools"
		samtools --version
		time samtools view -@ $THR -bS $TMP/tmp.$FROM.$TO/$FLNM.sam > $TMP/tmp.$FROM.$TO/$FLNM.bam
		rm $TMP/tmp.$FROM.$TO/$FLNM.sam
	else
		"Pair not found for $file"
	fi

	# updating iteration counter
	(( count++ ))
done < ~/tmp/R1.bwa.$FROM.$TO

echo "Date and time end of chunk pair $FROM to $TO out of 540 pairs"
date

echo "#== End of chunk for pairs $FROM to $TO ==#"
