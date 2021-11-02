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

# Run FastQC on 20 files at a time
cd $PATH_ORIG # change to raw fastqs dir
for i in `seq 0 20 1060` # 20 less than 1080!
do
	j=$((i+1)) # set upper boundary
	k=$((i+20)) # set lower boundary
	#echo "batch $j to $k" # just testing
	FLS=$(sed -n $j,${k}p ~/tmp/all_fastq_1080) #select 20 files
	echo "processing files $j to $k ..."
	time $FASTQC_PATH/fastqc -t 20 -out $OUT $FLS
	echo "The set of files $j to $k was succesfully analyzed with FastQC"	
	echo "The following files were analyzed:"
	echo $FLS
done

# This link on how to get a range of lines
## https://stackoverflow.com/questions/83329/how-can-i-extract-a-predetermined-range-of-lines-from-a-text-file-on-unix

