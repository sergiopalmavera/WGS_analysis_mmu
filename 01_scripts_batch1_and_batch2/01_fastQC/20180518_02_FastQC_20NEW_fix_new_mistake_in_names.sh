# I used "fastqc" instead of "fastq" to rename the reports

# Set the path to the results folder having the incorrect file names
PATH_RES=/projekte/I2-SOS-FERT/01_quality_control_fastqc/res

# Get the zip files that need renaming
ls -1 $PATH_RES | grep "001fastq-NEW" | grep ".zip" > ~/FBN_HOME/tmp/rename_zip.txt



# loop over the zip files
while read line
do
	NEW_NAME=$(echo $line | sed 's/001fastq/001.fastq/g')
	mv $PATH_RES/$line $PATH_RES/$NEW_NAME
done < ~/FBN_HOME/tmp/rename_zip.txt



# Get the html files that need renaming
ls -1 $PATH_RES | grep "N001fastq-NEW" | grep ".html" > ~/FBN_HOME/tmp/rename_html.txt

# loop over the html files
while read line
do	
	NEW_NAME=$(echo $line | sed 's/001fastq/001.fastq/g')
	mv $PATH_RES/$line $PATH_RES/$NEW_NAME
done < ~/FBN_HOME/tmp/rename_html.txt

