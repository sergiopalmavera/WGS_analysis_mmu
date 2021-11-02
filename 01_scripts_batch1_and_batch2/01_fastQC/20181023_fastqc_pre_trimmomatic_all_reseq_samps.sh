#!/bin/bash

# Define paths
ORIG=/projekte/I2-SOS-FERT/Original2
FASTQC=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC
OUTPUT=/projekte/I2-SOS-FERT/01_quality_control_fastqc2/res_not_concat

cd $ORIG
$FASTQC/fastqc -version
time $FASTQC/fastqc -t 20 -out $OUTPUT $(ls -1 | grep ".fastq.gz$")
