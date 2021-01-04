library(dplyr)
library(vroom)
library(stringr)
options(scipen=999)

# Define vcf file
tab <- "../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.table"

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
gt_tab <- vroom(tab) 

# chr X is turned into NAs
gt_tab$CHROM %>% unique()

# Fix chrX
gt_tab <- gt_tab %>% mutate(CHROM = ifelse(is.na(CHROM), "X", CHROM))

# Prepare col names
names(gt_tab) <- str_remove(names(gt_tab), ".GT") %>% str_remove("-L1")

# Match order of table and sample info
idx1 <- match(sample_info$sample_id, names(gt_tab))
tmp1 <- gt_tab[,idx1]
gt_tab <- gt_tab %>% 
  dplyr::select(CHROM,POS) %>% 
  bind_cols(tmp1) 

# Colnames match sample info?
identical(sample_info$sample_id,names(tmp1))

# Get sites with at least N non-missing samples
idx2 <- gt_tab %>% 
  dplyr::select(-CHROM,-POS) %>% 
  apply(1, function(x){

    #x=gt_tab %>% dplyr::select(-CHROM,-POS) %>% .[1,] %>% unlist()

    # count calls (non ./.) per group and check if N is >= 15 samples (except DU6)
    xx <- tapply(x[sample_info$Linie != "DU6"], 
           sample_info$Linie[sample_info$Linie != "DU6"], 
           function(gt) sum(gt != "./.") >= 15)
    # count calls (non ./.) in DU6 and check if N is >= 12 samples
    yy <- tapply(x[sample_info$Linie == "DU6"], 
                 sample_info$Linie[sample_info$Linie == "DU6"], 
                 function(gt) sum(gt != "./.") >= 12)
    
    # if all groups pass the min N per group (15-rest or 12-DU6)
    c(xx,yy) %>% all()
    
  })

sum(idx2)

# Define output file name
out_nm <- paste0("../output/keep_indels_DU6_12_REST_15", ".intervals")

# Get and export intervals for GATK
gt_tab %>% 
  dplyr::select(CHROM,POS) %>% 
  .[idx2,] %>% 
  mutate(tmp = paste0(CHROM,":",POS,"-",POS)) %>% 
  dplyr::select(tmp) %>% 
  write.table(out_nm, quote = F, col.names = F, row.names = F)
