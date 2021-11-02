#!/bin/bash

# Processing files in approx 1 day
## 1 day = 40 files (20 R1 , 20 R2)
## This corresponds to 20 pairs

# Define start and end of chunk
FROM=521
TO=540

# Define a counter variable to print out the number of iterations
count=$FROM

echo "Start of chunk for pairs $FROM to $TO"
echo "Defining paths"
ORIG=/projekte/I2-SOS-FERT/Original
ADAPTERS=/home/fb4/palma-vera/FBN_HOME/Tools/Trimmomatic-0.38/adapters
TRIMMOMATIC=/home/fb4/palma-vera/FBN_HOME/Tools/Trimmomatic-0.38
FASTQC=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC
OUTPUT=/projekte/I2-SOS-FERT/02_trimmed/results

echo "Defining R1 files for chunk"
ls -1 $ORIG | grep -i ".fastq.gz" | sed '/md5/d' | grep "R1" | sed -n $FROM,${TO}p > ~/tmp/R1.$FROM.$TO

echo "Date and time starting of chunk pairs $FROM to $TO out of 540 pairs"
date

echo "Creating temporary directory for chunk"
mkdir $OUTPUT/tmp.REDONE.$FROM.$TO
mkdir $OUTPUT/tmp.REDONE.$FROM.$TO/TRIMMO
mkdir $OUTPUT/tmp.REDONE.$FROM.$TO/FASTQC

while read file
do
	echo "#== Processing pair $count =="
	R1=$file
	R2=${file/R1/R2}
	echo $R1
	echo $R2

	echo "Looping over R1 and R2 pairs making first sure both files exist"
	if [ -e $ORIG/$R1 ] && [ -e $ORIG/$R2 ]
	then
		echo "Both files exist"
		echo "date and time sample pair started running through Trimmomatic and FastQC"
		date

		echo "Making sample directory name"
		TMP=$OUTPUT/tmp.REDONE.$FROM.$TO/${file%%_R*}
		mkdir -p $TMP

		echo "Making Trimmomatic and FastQC directories"
		mkdir -p $TMP/TRIMMO
		mkdir -p $TMP/FASTQC
		
		echo "Running Trimmomatic version"
		java -jar $TRIMMOMATIC/trimmomatic-0.38.jar -version
		time java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 $ORIG/$R1 $ORIG/$R2 $TMP/TRIMMO/${R1%%.fastq.gz}.OutputPaired.fastq.gz $TMP/TRIMMO/${R1%%.fastq.gz}.OutputUnpaired.fastq.gz $TMP/TRIMMO/${R2%%.fastq.gz}.OutputPaired.fastq.gz $TMP/TRIMMO/${R2%%.fastq.gz}.OutputUnpaired.fastq.gz ILLUMINACLIP:$ADAPTERS/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		echo "Running FastQC version"
		$FASTQC/fastqc -version
		THR=$(ls -1 $TMP/TRIMMO | wc -l)
		time $FASTQC/fastqc -t $THR -out $TMP/FASTQC $TMP/TRIMMO/*

		echo "Moving temporary results into temporary chunk directory"
		mv $TMP/TRIMMO/* $OUTPUT/tmp.REDONE.$FROM.$TO/TRIMMO
		mv $TMP/FASTQC/* $OUTPUT/tmp.REDONE.$FROM.$TO/FASTQC
		rm -rf $TMP

		echo "date and time pair finished running through Trimmomatic and FastQC"
		date
	else
		"Pair not found for $file"
	fi

	# updating iteration counter
	(( count++ ))
done < ~/tmp/R1.$FROM.$TO

echo "Date and time end of chunk pair $FROM to $TO out of 540 pairs"
date

echo "#== End of chunk for pairs $FROM to $TO ==#"
