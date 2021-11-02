FASTQ=/projekte/I2-SOS-FERT/02_trimmed2/results
OUT=/projekte/I2-SOS-FERT/02_trimmed2/GB_NR

FLS=$(ls -1 $FASTQ | grep "OutputPaired")

time for fastq in $FLS
do
	echo "## Calculating GB/read in  $fastq"
	time zcat $FASTQ/$fastq | paste - - - - | cut -f2 | tr -d '\n' | wc -c > $OUT/${fastq%.fastq.gz}.GB
	time zcat $FASTQ/$fastq | grep "@" | wc -l > $OUT/${fastq%.fastq.gz}.NR
done 
