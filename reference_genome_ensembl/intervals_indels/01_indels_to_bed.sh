#!/bin/bash

# Do not forget to habe bedops added in your PATH

bedops=/home/fb4/palma-vera/FBN_HOME/Tools/bedops/bin

echo "# Making bed file of deletions"
$bedops/vcf2bed --deletions < ../mus_musculus.vcf > ./mus_musculus_dels.bed
printf "\n"

echo "# Making bed file of insertions"
$bedops/vcf2bed --insertions < ../mus_musculus.vcf > ./mus_musculus_ins.bed
printf "\n"

echo "# Done"

