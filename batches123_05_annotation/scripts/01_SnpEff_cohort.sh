SnpEff=/home/fb4/palma-vera/FBN_HOME/Tools/SnpEff/snpEff
VARS=../../batches123_04_final_VCF/output
IN=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf
OUT=${IN/.vcf/.ann.vcf}

echo "# Stating SnpEff..."
java -Xmx4g -jar $SnpEff/snpEff.jar -csvStats ../output/${IN/.vcf/.summary.csv} GRCm38.86 $VARS/$IN > ../output/$OUT 
echo "#... SnpEff done"
