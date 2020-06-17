# Visualize the proportion of polymorphic, fixed-ref and fixed-alt
library(stringr)
library(ggplot2)
library(dplyr)

dat <- lapply(list.files("../output", pattern = "*.ALTfrq"), function(fl){
  
  pop <- fl %>% 
    str_remove(
      "cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered."
    ) %>% 
    str_remove(".ALTfrq")
  
  read.table(file.path("../output",fl)) %>% 
    mutate(Linie = pop)
}) %>% 
  bind_rows() 

head(dat)

dat <- dat %>% 
  mutate(state = NA) %>%
  mutate(state = ifelse(V1 > 0 & V1 < 1, "polymorphic",state)) %>% 
  mutate(state = ifelse(V1 == 0, "fixed_REF",state)) %>% 
  mutate(state = ifelse(V1 == 1, "fixed_ALT",state))

dat[sample(1:nrow(dat), 10),] %>% arrange(state) # corroborate

dat_summ <- dat %>% 
  group_by(Linie, state) %>% 
  summarise(n = n()) 

dat_summ %>% ungroup() %>% group_by(Linie) %>% summarise(sum(n)) # corroborate all records appear (5126585)

# visualize

dat_summ <- dat_summ %>% 
  ungroup() %>% 
  mutate(Linie = factor(Linie, levels = c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU")))


dat_summ %>% 
  ggplot(aes(x = Linie, y = n, fill = state)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    ylab("SNP counts") +
    xlab(NULL) +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle("Allele Frequency State Observed SNPs per Population") +
    ggsave("../figures/allele_frq_state.png")














