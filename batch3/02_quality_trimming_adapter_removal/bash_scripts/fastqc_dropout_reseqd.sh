#!/bin/bash

# Absolute paths (modify accordingly)
fastqc=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC

# Relative paths (leave as it is)
input=../output
output=../output_fastqc

# Define input files list

# Print program version
$fastqc/fastqc -version

# Run fastqc
time $fastqc/fastqc -t 50 -out $output $input/I34772-L1_S19_L004_R1_001.corrected.fastq $input/I34772-L1_S19_L004_R2_001.corrected.fastq

