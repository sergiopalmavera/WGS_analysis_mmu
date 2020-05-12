library(dplyr)
library(vroom)
library(stringr)

# Define vcf file
vcf <- "cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.table"

# sample info batch1-2 (60 samples, ~20x)
sample_info1 <- read.csv("../../sample_info/sample_info_batch1_batch2.csv", stringsAsFactors = F)

# sample info last batch (90 samples, ~5x avg-cvg)
sample_info2 <- read.csv("../../sample_info/sample_info_batch3/sample_info.csv", stringsAsFactors = F)

# Combine both data sets
ss1 <- sample_info1 %>% 
  dplyr::select(Linie, sample_id) %>% 
  mutate(target_cvg = "30x")

ss2 <- sample_info2 %>% dplyr::select(Linie, name) %>% 
  mutate(name = str_remove(name,"-S1")) %>% 
  dplyr::rename(sample_id = name) %>% 
  mutate(target_cvg = "5x",
         Linie = ifelse(Linie == "HLB", "DUhLB", Linie)) 

sample_info <- bind_rows(ss1,ss2) %>% 
  mutate(Linie = factor(Linie, c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU"))) %>% 
  arrange(Linie, sample_id) %>% 
  mutate(Linie = as.character(Linie))

# Load genotype table
gt_tab <- vroom(file.path("../output", vcf))

# Prepare col names
names(gt_tab) <- str_remove(names(gt_tab), ".GT") %>% str_remove("-L1")

# Match order of table and sample info
idx1 <- match(sample_info$sample_id, names(gt_tab))
tmp1 <- gt_tab[,idx1]
gt_tab <- gt_tab %>% 
  dplyr::select(CHROM,POS) %>% 
  bind_cols(tmp1) 

# Create index. TRUE for SNPs in which missingness is below 50%
idx2 <- gt_tab %>% 
  dplyr::select(-CHROM,-POS) %>% 
  apply(1, function(x){
    tapply(x, sample_info$Linie, function(gt) mean(gt != "./.") >= 0.9) %>% any()
  })

# Get and export intervals for GATK
gt_tab %>% 
  dplyr::select(CHROM,POS) %>% 
  .[idx2,] %>% 
  mutate(tmp = paste0(CHROM,":",POS,"-",POS)) %>% 
  dplyr::select(tmp) %>% 
  write.table("../output/keep_snps.intervals", quote = F, col.names = F, row.names = F)
