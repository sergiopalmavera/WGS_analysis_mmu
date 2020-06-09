library(dplyr)
library(vroom)

# Old VCF (10 animals per line, total 60 animals)
snps_vcf1 <- vroom(
  # absolute path
  "/projekte/I2-SOS-FERT/15_GenotypeRefinement/truth_sens_85pct/by_line/GroupVCFs_merged_minDPmaxDPminGQ/minDP5minGQ20maxDP/final_vcf_snp_pos.txt", 
  col_names = F
  )

# New VCF (includes 10 animals per line, plus 15 new animals per line. Total 150 animals)
snps_vcf2 <- vroom("batches123_04_FinalVCF/output/final_vcf_snp_pos.txt", col_names = F)


intx <- inner_join(snps_vcf1, snps_vcf2, by = c("X1","X2"))

dim(intx) #3508724       2
