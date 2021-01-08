GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

vcf=mgp.v5.merged.indels.dbSNP142.normed.vcf.gz
vcf_pass=${vcf/.vcf.gz/.PASS.vcf.gz}

$GATK/gatk IndexFeatureFile -F $vcf

$GATK/gatk SelectVariants \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V $vcf \
	-O $vcf_pass \
	-select-type INDEL \
	--exclude-filtered

echo "## Number of sites in new PASSonly vcf (gatk SelectVariants)"
grep -v "^#" $vcf_pass | wc -l


echo "## Corroborate only PASS sites in FILTER"
echo $vcf_pass
grep -v "^#" $vcf_pass | awk '{print $7}' | sort | uniq

