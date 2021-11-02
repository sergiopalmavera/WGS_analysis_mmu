# this script compares the snps discovered in the original cohort (60 genomes) vs the current vcf with 150 genomes.

library(dplyr)

snps_60_genomes <- read.table(
  "/projekte/I2-SOS-FERT/15_GenotypeRefinement/truth_sens_85pct/by_line/GroupVCFs_merged_minDPmaxDPminGQ/minDP5minGQ20maxDP/final_vcf_snp_pos.txt", 
  header = F, 
  sep = ":", 
  stringsAsFactors = F
  )

head(snps_60_genomes)

snps_150_genomes <- read.table(
  "../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.SNPids",
  header = F,
  stringsAsFactors = F
  )

head(snps_150_genomes)

intx <- inner_join(
  snps_60_genomes,
  snps_150_genomes,
  by = c("V1","V2")
)

dim(intx)

nrow(intx)/nrow(snps_60_genomes)

nrow(intx)/nrow(snps_150_genomes)
