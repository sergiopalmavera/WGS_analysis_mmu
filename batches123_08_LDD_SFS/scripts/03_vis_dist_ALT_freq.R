library(dplyr)
library(ggplot2)
library(stringr)

dat <- lapply(list.files("../output", pattern = "*.ALTfrq"), function(fl){
  
  pop <- fl %>% 
    str_remove(
      "cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.allrecords."
      ) %>% 
    str_remove(".ALTfrq")
  
  read.table(file.path("../output",fl)) %>% 
    mutate(Linie = pop)
}) %>% 
  bind_rows() 

dat <- dat %>% mutate(Linie = factor(Linie, c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU")))

dat %>% 
  ggplot(aes(x = V1)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Linie, ncol = 2) +
    theme_bw() +
    xlab("Frequency") +
    ggtitle("Alternative Allele Frequency") +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    ggsave("../figures/hist_ALT_frq.png")
