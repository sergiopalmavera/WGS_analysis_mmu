#!/bin/bash

# Processing files in approx 1 day
## 1 day = 40 files (20 R1 , 20 R2)
## This corresponds to 20 pairs

# Define start and end of chunk
FROM=1
TO=20

# Define paths
ORIG=/projekte/I2-SOS-FERT/Original
ADAPTERS=/home/fb4/palma-vera/FBN_HOME/Tools/Trimmomatic-0.38/adapters
TRIMMOMATIC=/home/fb4/palma-vera/FBN_HOME/Tools/Trimmomatic-0.38
FASTQC=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC
OUTPUT=/projekte/I2-SOS-FERT/02_trimmed/results

# Define R1 files for chunk
ls -1 $ORIG | grep -i ".fastq.gz" | sed '/md5/d' | grep "R1" | sed -n $FROM,${TO}p > ~/tmp/R1.$FROM.$TO

while read file
do
	echo "date and time chunk started"
	date

	echo "Processing pair of files:"
	R1=$file
	R2=${file/R1/R2}
	echo $R1
	echo $R2

	echo "Looping over R1 and R2 pairs making first sure both files exist"
	if [ -e $ORIG/$R1 ] && [ -e $ORIG/$R2 ]
	then
		echo "Both files exist"

		echo "Making root directory name"
		TMP=${file%%_R*}

		echo "Making temporary directories"
		mkdir -p $OUTPUT/$TMP
		mkdir -p $OUTPUT/$TMP/TRIMMO
		mkdir -p $OUTPUT/$TMP/FASTQC
		
		echo "Running Trimmomatic"
		time java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 $ORIG/$R1 $ORIG/$R2 $OUTPUT/$TMP/TRIMMO/${R1%%.fastq.gz}.OutputPaired.fastq.gz $OUTPUT/$TMP/TRIMMO/${R1%%.fastq.gz}.OutputUnpaired.fastq.gz $OUTPUT/$TMP/TRIMMO/${R2%%.fastq.gz}.OutputPaired.fastq.gz $OUTPUT/$TMP/TRIMMO/${R2%%.fastq.gz}.OutputUnpaired.fastq.gz ILLUMINACLIP:$ADAPTERS/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		echo "Running FastQC"
		THR=$(ls -1 $OUTPUT/$TMP/TRIMMO | wc -l)
		time $FASTQC/fastqc -t $THR -out $OUTPUT/$TMP/FASTQC $OUTPUT/$TMP/TRIMMO/*

		echo "Moving temporary results into permanent directory"
		mv $OUTPUT/$TMP/TRIMMO/* $OUTPUT/TRIMMO
		mv $OUTPUT/$TMP/FASTQC/* $OUTPUT/FASTQC
		rm -rf $OUTPUT/$TMP
	else
		"Pair not found for $file"
	fi
	echo "Date and time of end"
	date
done < ~/tmp/R1.$FROM.$TO
