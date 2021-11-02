#!/bin/bash

chr=$1

sed 's/chr=/chr='$1'/' GenotypeGVCFs_template.sh > GenotypeGVCFs_chr${chr}.sh

chmod +x GenotypeGVCFs_chr${chr}.sh

nohup GenotypeGVCFs_chr${chr}.sh &> GenotypeGVCFs_chr${chr}.out &

sleep 20

ps -xjf | head -1 > tmp1
ps -xjf | grep GenotypeGVCFs_chr${chr}.sh | grep bash > tmp2

cat tmp1 tmp2 > GenotypeGVCFs_chr${chr}.pid

rm tmp1 tmp2

