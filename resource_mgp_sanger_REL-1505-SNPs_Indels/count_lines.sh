echo "## Total lines in VCF before any filtering (incl SNPS, indels, PASS and not PASS)" 
echo "## File mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz | grep -v "^#" | wc -l

echo "## Number of sites in new PASSonly vcf (gatk SelectVariants)"
echo "## File mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf" 
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf | wc -l

echo "## Number of sites in previous PASSonly vcf (VCFtools)"
echo "## File:mgp.v5.merged.snps_all.dbSNP142_PASS.vcf"
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS.vcf | wc -l

echo "## Corroborate only PASS sites in FILTER"

echo "mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf"
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS_final.vcf | awk '{print $7}' | sort | uniq

echo "mgp.v5.merged.snps_all.dbSNP142_PASS.vcf"
grep -v "^#" mgp.v5.merged.snps_all.dbSNP142_PASS.vcf | awk '{print $7}' | sort | uniq
