library(here)
library(dplyr)

sample_info <- read.csv(here("sample_info/sample_info_batch1_batch2.csv")) %>% 
  dplyr::select(Linie, sample_id) %>% 
  mutate(target_cvg = "30x") %>% 
  bind_rows(
    read.csv(here("sample_info/sample_info_batch3/sample_info.csv")) %>% 
      dplyr::select(Linie,name) %>% 
      dplyr::rename(sample_id = name)
  ) %>% 
  mutate(Linie = ifelse(Linie == "HLB", "DUhLB",Linie))



final_vcf_metrics_detail <- read.delim(here("batches123_04_FinalVCF/metrics/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.table.variant_calling_detail_metrics"), stringsAsFactors = F, comment.char = "#", dec = ",") %>% 
  mutate(SAMPLE_ALIAS2 = str_remove(SAMPLE_ALIAS, "-L1")) %>% 
  left_join(sample_info, by = c("SAMPLE_ALIAS2" = "sample_id")) %>% 
  #mutate(Linie = factor(Linie, levels = c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU"))) %>% 
  arrange(Linie, target_cvg, SAMPLE_ALIAS2)



final_vcf_metrics_detail$Linie %>% unique()


fx <- function(l){
  final_vcf_metrics_detail %>% 
    filter(Linie == l) %>% 
    dplyr::select(SAMPLE_ALIAS) %>% 
    write.table(paste0("sample_info/vcf_samples_",l), quote = F, col.names = F, row.names = F)  
}

fx("DUK")
fx("DUC")
fx("DU6")
fx("DU6P")
fx("DUhLB")
fx("FZTDU")
