
# Define path to original files
PATH_ORIG=/projekte/I2-SOS-FERT/Original

# Define the files that need to be processed and export them:
ls $PATH_ORIG | grep fastq-NEW.gz > ~/tmp/NEW_FASTQs.txt # generic name, details in .out file

## Set path to FastQC
FASTQC_PATH=/home/fb4/palma-vera/Tools/FastQC

## Set path for output
OUT=/projekte/I2-SOS-FERT/01_quality_control_fastqc/res

## make a file to store results for renaming
RENAME=/projekte/I2-SOS-FERT/01_quality_control_fastqc/tmp/rename
rm -r $RENAME # in case it was already made while testing
mkdir $RENAME

## Add path to files 
while read line
do
	echo $PATH_ORIG/$line
done < ~/tmp/NEW_FASTQs.txt > ~/tmp/NEW_FASTQs_fullpaths.txt

## put paths and files in a variable
FLS=$(cat ~/tmp/NEW_FASTQs_fullpaths.txt)

## Run FastQC in all 20 samples in parallel
time $FASTQC_PATH/fastqc -t 20 -out $RENAME $FLS

## print files 
echo "Files processed:"
cat ~/tmp/NEW_FASTQs.txt

