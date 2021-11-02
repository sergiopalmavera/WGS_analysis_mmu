#!/bin/bash

FASTQ=../../02_quality_trimming_adapter_removal/output
OUT=../output_n_bp
TMP=$OUT/TMP

R1=I34772-L1_S19_L004_R1_001.corrected.fastq
R2=I34772-L1_S19_L004_R2_001.corrected.fastq

time cat $FASTQ/$R1 | paste - - - - | cut -f2 | tr -d '\n' | wc -c > $OUT/${R1%.fastq}.NBP
time cat $FASTQ/$R1 | grep "@" | wc -l > $OUT/${R1%.fastq}.NR

time cat $FASTQ/$R2 | paste - - - - | cut -f2 | tr -d '\n' | wc -c > $OUT/${R2%.fastq}.NBP
time cat $FASTQ/$R2 | grep "@" | wc -l > $OUT/${R2%.fastq}.NR

