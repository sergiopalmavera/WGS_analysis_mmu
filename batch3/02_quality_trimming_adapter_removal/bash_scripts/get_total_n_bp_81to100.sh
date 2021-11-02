#!/bin/bash

FROM=81
TO=100

FASTQ=../../02_quality_trimming_adapter_removal/output
OUT=../output_n_bp
TMP=$OUT/TMP

FLS=$(ls -1 $FASTQ | grep ".fastq.gz" | sed -n $FROM,${TO}p)
echo "# Processing files"
for fastq in $FLS; do echo $fastq; done

time for fastq in $FLS
do
	echo "## Calculating GB/read in  $fastq"
	time zcat $FASTQ/$fastq | paste - - - - | cut -f2 | tr -d '\n' | wc -c > $OUT/${fastq%.fastq.gz}.NBP
	time zcat $FASTQ/$fastq | grep "@" | wc -l > $OUT/${fastq%.fastq.gz}.NR
done 
