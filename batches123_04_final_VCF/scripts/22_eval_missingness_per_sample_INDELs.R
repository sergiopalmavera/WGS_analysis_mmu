library(vroom)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)

# name of data file
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

# Prepare col names
names(gt_tab) <- str_remove(names(gt_tab), ".GT") %>% str_remove("-L1")
dim(gt_tab)
 
# Match order of table and sample info
idx1 <- match(sample_info$sample_id, names(gt_tab))
tmp1 <- gt_tab[,idx1]
gt_tab <- gt_tab %>% 
   dplyr::select(CHROM,POS) %>% 
   bind_cols(tmp1) 
 
# columns in the right order?
identical(sample_info$sample_id,names(tmp1))
 
# make a function to extract number of missing samples per pop at each SNP
get_n_not_miss <- function(l){
   i_col <- names(gt_tab) %in% sample_info$sample_id[sample_info$Linie == l]
   gt_tab[,i_col] %>% apply(1,function(row_x) sum(row_x != "./."))
}
 
# Put count missings per populations counts into one data frame
not_miss_cts <- data.frame(
   DUK = get_n_not_miss("DUK"),
   DUC = get_n_not_miss("DUC"),
   DU6 = get_n_not_miss("DU6"),
   DU6P = get_n_not_miss("DU6P"),
   DUhLB = get_n_not_miss("DUhLB"),
   FZTDU = get_n_not_miss("FZTDU")
 ) 
 
# define min number of called-samples in pop 
minNs <- 0:25
 
# Count number of SNPs per population with at least minN calls
n_indels <- data.frame(
   min_called_samples = minNs,
   DUK = vector(mode = "numeric", length = length(minNs)),
   DUC = vector(mode = "numeric", length = length(minNs)),
   DU6 = vector(mode = "numeric", length = length(minNs)),
   DU6P = vector(mode = "numeric", length = length(minNs)),
   DUhLB = vector(mode = "numeric", length = length(minNs)),
   FZTDU = vector(mode = "numeric", length = length(minNs))
 )
 
for(i in 1:length(minNs)){
   
   minN <- minNs[i]
   
   res <- sapply(not_miss_cts, function(col_i) sum(col_i >= minN) ) 
   
   n_indels[i,2:ncol(n_indels)] <- res
 }
 
# convert to fraction
fraction_indels <- lapply(n_indels[2:ncol(n_indels)], function(x) x/nrow(gt_tab)) %>% 
   bind_cols() %>% 
   mutate(min_called_samples = minNs) %>% 
   dplyr::select(min_called_samples, everything())
 
 
write.csv(n_indels, "../figures_tables/n_indels_by_min_n_called_per_line.csv")
 
write.csv(fraction_indels, "../figures_tables/fraction_indels_by_min_n_called_per_line.csv")
