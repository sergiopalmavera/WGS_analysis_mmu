RES=/projekte/I2-SOS-FERT/02_trimmed/results

ls -1 $RES | grep -i "tmp" > ~/tmp/tmp.files

while read name
do
	echo "# Temp directory: $name ---------------------"
	echo
	echo "Moving Trimmomatic files"
	mv $RES/$name/TRIMMO/* $RES/TRIMMO
	echo
	echo "Moving FastQC files"
	mv $RES/$name/FASTQC/* $RES/FASTQC
	echo
	echo "Removing temporary folder"
	rm -r $RES/$name
done < ~/tmp/tmp.files
