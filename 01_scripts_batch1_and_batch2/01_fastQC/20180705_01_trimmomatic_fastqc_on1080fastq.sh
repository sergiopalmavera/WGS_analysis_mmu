##== Description ==
# Processing all raw fastqs with trimmomatic, followed by FastQC.
# Raw reads had the problem that the FastQC modules "Kmer content" and "overrepresented sequences" failed or gave warnings. 
# I did a test to remove this issues with trimmomatic.
# For that I used the adapter file provided by trimmomatic (TruSeq3-PE-2.fa). "overrepresented sequences" was fixed, but "kmer content" was not. However, "kmer content" graph changed.
# I posted the question on Biostars and it is apparently a minor problem, and it should be good enough to move on.

echo "##== Define paths ================"
ORIG=/projekte/I2-SOS-FERT/Original
ADAPTERS=/home/fb4/palma-vera/FBN_HOME/Tools/Trimmomatic-0.38/adapters
TRIMMOMATIC=/home/fb4/palma-vera/FBN_HOME/Tools/Trimmomatic-0.38
FASTQC=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC
OUTPUT=/projekte/I2-SOS-FERT/02_trimmed/results

echo "##== Making list of R1 file names =========================="
ls -1 $ORIG | grep -i ".fastq.gz" | sed '/md5/d' | grep "R1" > ~/tmp/R1

echo "##== Looping over R1 and R2 pairs =========================="
while read file
do 
	echo "##== Capturing Staring Date =============================="
	date

	echo "##== Processing pairs of files =================="
	R1=$file
	R2=${file/R1/R2}
	echo $R1
	echo $R2
	
	echo "##== Looping over R1 and R2 pairs making sure both files exist  ==========="
	if [ -e $ORIG/$R1 ] && [ -e $ORIG/$R2 ]
	then
		echo "##== Both files exist! ==========="
		echo
		echo "##== Making temporary directories ====="
		mkdir -p $OUTPUT/tmp
		mkdir -p $OUTPUT/tmp/TRIMMO
		mkdir -p $OUTPUT/tmp/FASTQC
		
		echo "##== Running Trimmomatic =========="
		time java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 $ORIG/$R1 $ORIG/$R2 $OUTPUT/tmp/TRIMMO/${R1%%.fastq.gz}.OutputPaired.fastq.gz $OUTPUT/tmp/TRIMMO/${R1%%.fastq.gz}.OutputUnpaired.fastq.gz $OUTPUT/tmp/TRIMMO/${R2%%.fastq.gz}.OutputPaired.fastq.gz $OUTPUT/tmp/TRIMMO/${R2%%.fastq.gz}.OutputUnpaired.fastq.gz ILLUMINACLIP:$ADAPTERS/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		echo "##== Running FastQC ==============="
		THR=$(ls -1 $OUTPUT/tmp/TRIMMO | wc -l)	
		time $FASTQC/fastqc -t $THR -out $OUTPUT/tmp/FASTQC $OUTPUT/tmp/TRIMMO/*

		echo "##== Moving results from tmp dir into permanent directory ========="
		mv $OUTPUT/tmp/TRIMMO/* $OUTPUT/TRIMMO	
		mv $OUTPUT/tmp/FASTQC/* $OUTPUT/FASTQC	
		rm -rf $OUTPUT/tmp
	else
		echo "Pair R2 file missing for $file"
		echo "Find out what is happening with $file pair"
	fi
	
	echo "##== Capturing Ending Time Stamp  ================================"
	date
done < ~/tmp/R1

##== Result ===============
