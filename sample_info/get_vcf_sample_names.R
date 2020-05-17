sample_info <- read.csv(here("sample_info/sample_info_batch1_batch2.csv")) %>% 
  dplyr::select(Linie, sample_id) %>% 
  mutate(target_cvg = "30x") %>% 
  bind_rows(
    read.csv(here("sample_info/sample_info_batch3/sample_info.csv")) %>% 
      dplyr::select(Linie,name) %>% 
      dplyr::rename(sample_id = name)
  ) %>% 
  mutate(Linie = ifelse(Linie == "HLB", "DUhLB",Linie))


sample_info %>% head()

sample_info$Linie %>% unique()


fx <- function(l){
  sample_info %>% filter(Linie == l) %>% dplyr::select(sample_id) %>% 
    write.table(paste0("sample_info/vcf_samples_",l), quote = F, col.names = F, row.names = F)  
}

fx("DUK")
fx("DUC")
fx("DU6")
fx("DU6P")
fx("DUhLB")
fx("FZTDU")
