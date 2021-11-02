#!/bin/bash

# QC raw reads of drop-out sample that was resequenced

# Absolute paths (modify accordingly)
input=/projekte/I2-SOS-FERT/Original3
fastqc=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC

# Relative paths (leave as it is)
output=../output


# Print program version
$fastqc/fastqc -version

# Run fastqc
time $fastqc/fastqc -t 50 -out $output $input/I34772-L1_S19_L004_R1_001.fastq.gz $input/I34772-L1_S19_L004_R2_001.fastq.gz 

