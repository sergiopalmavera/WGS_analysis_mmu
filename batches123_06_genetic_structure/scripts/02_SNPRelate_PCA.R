# Checking example before running SNPRelate on my data
## According to https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate
## To clarify axes labeling: https://support.bioconductor.org/p/119389/

# Load packages
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(dplyr)

# convert VCF to GDS
vcf <- "../data/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ldpruned.vcf"

seqVCF2GDS(vcf, "../data/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ldpruned.gds")

# Import GDS data
genofile <- seqOpen("../data/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ldpruned.gds")

# PCA
pca <- snpgdsPCA(genofile)

# Export table of PC percentage
data.frame(PC=1:7, varprop = pca$varprop[1:7]*100) %>% write.csv("../data/pc_pct.csv")

# Export table of eigenscores
tab <- data.frame(sample.id = pca$sample.id,
		  PC1 = pca$eigenvect[,1],    # the first eigenvector
		  PC2 = pca$eigenvect[,2],    # the second eigenvector
		  PC3 = pca$eigenvect[,3],    # the third eigenvector 
		  PC4 = pca$eigenvect[,4],    # the fourth eigenvector 
		  PC5 = pca$eigenvect[,5],
		  PC6 = pca$eigenvect[,6],
		  PC7 = pca$eigenvect[,7],
		  stringsAsFactors = FALSE)

write.csv(tab, "../data/pc_eigenscores.csv")

sessionInfo()
