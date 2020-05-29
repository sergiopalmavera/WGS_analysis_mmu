# All sites in cohort vcf were used.

library(ggplot2)
library(stringr)

dat <- lapply(list.files("../output", pattern = "*.plink.maf"), function(fl){
  
  pop <- fl %>% 
    str_remove(
      "cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.allrecords."
    ) %>% 
    str_remove(".plink.maf")
  
  read.table(file.path("../output",fl), header = T) %>% 
    mutate(Linie = pop)
}) %>% 
  bind_rows() 

dat <- dat %>% mutate(Linie = factor(Linie, c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU")))

head(dat)

dat %>% 
  # remove MAF=0
  filter(MAF > 0) %>% 
  ggplot(aes(x = MAF)) +
  geom_histogram(bins = 50) +
  facet_wrap(~Linie, ncol = 2) +
  theme_bw() +
  xlab("Frequency") +
  ggtitle("Minor Allele Frequency (excl MAF=0)") +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggsave("../figures/hist_MAF_excl_maf0.png")

dat %>% 
  ggplot(aes(x = MAF)) +
  geom_histogram(bins = 50) +
  facet_wrap(~Linie, ncol = 2) +
  theme_bw() +
  xlab("Frequency") +
  ggtitle("Minor Allele Frequency (incl MAF=0)") +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggsave("../figures/hist_MAF_incl_maf0.png")

