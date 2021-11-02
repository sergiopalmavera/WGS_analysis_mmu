# This script trimms reads according to quality, detects and removes adatpters.
# https://github.com/OpenGene/fastp
# This software solves the problem with polyG tails inherent to NovaSeq and not found in HiSeq reads.
# check this post of mine to learn more about this issue: https://www.biostars.org/p/359760/#359822

# Absolute paths (modify accordingly)
input=/projekte/I2-SOS-FERT/Original3
fastp=~/FBN_HOME/Tools/fastp

# Relative paths (leave as it is)
output=../output

# Define files to process, replacing read pair with RX as placeholder ==> pairs 
fls=$(ls -1 $input/*.fastq.gz | sed 's/R[1-2]/RX/' | sort | uniq)

# Loop over each distinct sample and run fastp on the corresponding read pair
for fl in $fls
do
	echo "# Starting with pair $fl"

	# Define read file names
	R1=$(echo $fl | sed 's/RX/R1/')
	R2=$(echo $fl | sed 's/RX/R2/')
	
	# Define output file names
	nm_html=$(basename $fl | sed 's/RX/[R1,R2]/; s/.fastq.gz/.html/')
	nm_json=$(basename $fl | sed 's/RX/[R1,R2]/; s/.fastq.gz/.json/')
	nm_out_R1=$(basename $R1 | sed 's/.fastq/.corrected.fastq/')
	nm_out_R2=$(basename $R2 | sed 's/.fastq/.corrected.fastq/')

	# Runf fastp 
	time $fastp/fastp -h $output/$nm_html -j $output/$nm_json -i $R1 -I $R2 -o $output/$nm_out_R1 -O $output/$nm_out_R2
	
	echo "done"
	printf "\n\n"
done




