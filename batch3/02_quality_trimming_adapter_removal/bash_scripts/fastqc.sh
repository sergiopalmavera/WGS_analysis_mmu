#!/bin/bash

# Absolute paths (modify accordingly)
fastqc=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC

# Relative paths (leave as it is)
input=../output
output=../output_fastqc

# Define input files list
fls=$(ls -1 $input/*.fastq.gz)

# Print program version
$fastqc/fastqc -version

# Run fastqc
time $fastqc/fastqc -t 50 -out $output $fls

