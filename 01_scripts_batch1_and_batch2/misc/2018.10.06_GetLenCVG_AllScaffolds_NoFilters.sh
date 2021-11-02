IN=/projekte/I2-SOS-FERT/03_alignments_raw/DepthAlongGenome_NoFilters/results
OUT=/projekte/I2-SOS-FERT/03_alignments_raw/DepthAlongGenome_NoFilters/res_for_R

cd $IN

time for FL in * 
do
	echo "# FILE: $FL"
	echo "## LEN CVG ALL SCAFFOLDS $FL"
	LEN=$(wc -l $FL | awk '{print $1}')
	echo $LEN
	#LEN=$(tail $FL |wc -l)
	for thr in `seq 1 4 49`
	do
		echo "### THRESHOLD $thr"
		CVG=$(awk -v x="$thr" '{if($3>=x) print}' $FL | wc -l)
		echo $CVG
		#CVG=$(tail $FL | awk -v x="$thr" '{if($3>=x) print}' | wc -l)
		RES=$(bc -l <<< "scale=2; $CVG/$LEN")
		echo $thr $RES $LEN
	done > $OUT/${FL/.cvg/_LenCVG_AllScaff}.txt
done
