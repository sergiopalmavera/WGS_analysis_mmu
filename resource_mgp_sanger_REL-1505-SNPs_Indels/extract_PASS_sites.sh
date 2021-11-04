GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=/projekte/I2-SOS-FERT/reference_genome_ensembl

$GATK/gatk SelectVariants \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
	-O mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf \
	-select-type SNP \
	--exclude-filtered

echo "## Number of sites in new PASSonly vcf (gatk SelectVariants)"
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf | wc -l

echo "## Number of sites in previous PASSonly vcf (VCFtools)"
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS.vcf | wc -l

echo "## Corroborate only PASS sites in FILTER"
echo "mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf"
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf | awk '{print $7}' | sort | uniq

echo "mgp.v5.merged.snps_all.dbSNP142_PASS.vcf"
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS.vcf | awk '{print $7}' | sort | uniq
