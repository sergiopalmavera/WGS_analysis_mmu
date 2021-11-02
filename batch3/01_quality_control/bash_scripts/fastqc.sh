#!/bin/bash

# Absolute paths (modify accordingly)
input=/projekte/I2-SOS-FERT/Original3
fastqc=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC

# Relative paths (leave as it is)
output=../output

# Define input files list
fls=$(ls -1 $input | grep ".fastq.gz$" | for i in $(cat); do echo "$input/$i";done)

# Print program version
$fastqc/fastqc -version

# Run fastqc
time $fastqc/fastqc -t 50 -out $output $fls

