# Files are located in the folder
RENAME=/projekte/I2-SOS-FERT/01_quality_control_fastqc/tmp/rename

# They need to be renamed and moved to the output folder
OUT=/projekte/I2-SOS-FERT/01_quality_control_fastqc/res


# Loop over each file, renaming and moving it to the final folder
for file in $RENAME/*
do
	# remove unwanted part of name
	NEW=${file:59}
	NEW=$(echo "${NEW//fastq-}" ) 

	# move file into final folder
	mv $file $OUT/$NEW 
done


