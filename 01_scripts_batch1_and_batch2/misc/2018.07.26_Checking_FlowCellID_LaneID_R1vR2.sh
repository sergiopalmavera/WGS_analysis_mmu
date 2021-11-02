# Define paths
ORIG=/projekte/I2-SOS-FERT/Original

# Get first raw fastqR1
while read R1
do
	SAMP=$(echo $R1 | cut -d '_' -f1) 
	R2=${R1/R1/R2}
	echo "Sample: $SAMP"
	# get flow cell info
	FLOWCELLR1=$(zcat $ORIG/$R1 | grep '@' | cut -d ':' -f3 | sort | uniq)
	FLOWCELLR2=$(zcat $ORIG/$R2 | grep '@' | cut -d ':' -f3 | sort | uniq)
 	# get lane info
	LANER1=$(zcat $ORIG/$R1 | grep '@' | cut -d ':' -f4 | sort | uniq)
	LANER2=$(zcat $ORIG/$R2 | grep '@' | cut -d ':' -f4 | sort | uniq)
	# check flow cells are equal
	if [ $FLOWCELLR1 = $FLOWCELLR2 ]
	then
		echo "Flow cells equal"
	else
		echo "Flow cells not equal"
	fi
	# check lanes are equal
	if [ $LANER1 = $LANER2 ]
	then
		echo "Lanes equal"
	else
		echo "Lanes not equal"
	fi
	# check both 
	if [  $FLOWCELLR1 = $FLOWCELLR2 ] && [ $LANER1 = $LANER2 ]
	then
		echo "flow cells and lanes are same"
	else
		echo "mismatches in lane or flowcell"
	fi
done <~/tmp/R1all540.txt
