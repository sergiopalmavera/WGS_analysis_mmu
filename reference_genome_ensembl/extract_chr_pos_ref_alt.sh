#/bin/bash
vcf=mus_musculus.vcf
bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
$bcftools/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf > ${vcf/.vcf/.tab}
