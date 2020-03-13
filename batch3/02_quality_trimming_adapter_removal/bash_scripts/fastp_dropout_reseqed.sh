# This script trimms reads according to quality, detects and removes adatpters.
# https://github.com/OpenGene/fastp
# This software solves the problem with polyG tails inherent to NovaSeq and not found in HiSeq reads.
# check this post of mine to learn more about this issue: https://www.biostars.org/p/359760/#359822

# Absolute paths (modify accordingly)
input=/projekte/I2-SOS-FERT/Original3
fastp=~/FBN_HOME/Tools/fastp

# Relative paths (leave as it is)
output=../output

time $fastp/fastp -h $output/I34772-L1_S19_L004_[R1,R2]_001.html -j $output/I34772-L1_S19_L004_[R1,R2]_001.json -i $input/I34772-L1_S19_L004_R1_001.fastq.gz -I $input/I34772-L1_S19_L004_R2_001.fastq.gz -o $output/I34772-L1_S19_L004_R1_001.corrected.fastq -O $output/I34772-L1_S19_L004_R2_001.corrected.fastq
