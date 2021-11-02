#!/bin/bash

# Processing files in approx 1 day
## 1 day = 40 files (20 R1 , 20 R2)
## This corresponds to 20 pairs

# Define start and end of chunk
FROM=61
TO=80

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

		echo "Making root directory name"
		TMP=${file%%_R*}

		echo "Making temporary directories"
		mkdir -p $OUTPUT/$TMP
		mkdir -p $OUTPUT/$TMP/TRIMMO
		mkdir -p $OUTPUT/$TMP/FASTQC
		
		echo "Running Trimmomatic version"
		java -jar $TRIMMOMATIC/trimmomatic-0.38.jar -version
		time java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 $ORIG/$R1 $ORIG/$R2 $OUTPUT/$TMP/TRIMMO/${R1%%.fastq.gz}.OutputPaired.fastq.gz $OUTPUT/$TMP/TRIMMO/${R1%%.fastq.gz}.OutputUnpaired.fastq.gz $OUTPUT/$TMP/TRIMMO/${R2%%.fastq.gz}.OutputPaired.fastq.gz $OUTPUT/$TMP/TRIMMO/${R2%%.fastq.gz}.OutputUnpaired.fastq.gz ILLUMINACLIP:$ADAPTERS/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		echo "Running FastQC version"
		$FASTQC/fastqc -version
		THR=$(ls -1 $OUTPUT/$TMP/TRIMMO | wc -l)
		time $FASTQC/fastqc -t $THR -out $OUTPUT/$TMP/FASTQC $OUTPUT/$TMP/TRIMMO/*

		echo "Moving temporary results into permanent directory"
		mv $OUTPUT/$TMP/TRIMMO/* $OUTPUT/TRIMMO
		mv $OUTPUT/$TMP/FASTQC/* $OUTPUT/FASTQC
		rm -rf $OUTPUT/$TMP

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
