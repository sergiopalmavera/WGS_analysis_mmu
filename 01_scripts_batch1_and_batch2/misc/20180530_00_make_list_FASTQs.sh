# Intro
## Better run all raw FASTQs in one job
## Previous results were removed and corresponding directories cleared

# Set path to raw FASTQs
PATH_ORIG=/projekte/I2-SOS-FERT/Original

# Set path to FastQC
FASTQC_PATH=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC

# Set path for output
OUT=/projekte/I2-SOS-FERT/01_quality_control_fastqc/res

# Make a list with  all FASTQs
echo "The number of FASTQ files is: "
ls -1 $PATH_ORIG | grep "fastq" | sed '/md5/d' | wc -l  
ls -1 $PATH_ORIG | grep "fastq" | sed '/md5/d' > ~/tmp/all_fastq_1080
