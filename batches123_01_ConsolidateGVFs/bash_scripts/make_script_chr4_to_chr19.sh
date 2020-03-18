for i in $(seq 4 1 19)
do
	cp ConsolidateGVCFs_chr1.sh ./tmp/tmp${i} 
	sed "s/chr=1/chr=$i/" ./tmp/tmp${i} > ./ConsolidateGVCFs_chr${i}.sh
	nohup ./ConsolidateGVCFs_chr${i}.sh &> ./ConsolidateGVCFs_chr${i}.out &
	rm tmp/*
done
