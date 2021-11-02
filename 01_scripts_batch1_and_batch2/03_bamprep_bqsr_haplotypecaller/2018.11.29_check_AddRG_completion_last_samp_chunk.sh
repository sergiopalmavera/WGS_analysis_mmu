#!/bin/bash

FILE=$1

cd ~/FBN_HOME/Scripts_and_files/Jobs

LAST=$(sed -n 3p $FILE | awk '{print $10}')

LINE=$(grep -n $LAST $FILE | grep "Processing bam files for sample" | awk -F: '{print $1}')

Nfiles=$(tail --lines=+$LINE $FILE | grep "Number of bams for sample" | awk '{print $7}')

Ndone=$(tail --lines=+$LINE $FILE | grep done | wc -l)

if [ $Nfiles -eq $Ndone ]
then 
	echo "RG added; ready for merging"
	echo "Last sample in chunk: $LAST"
fi

