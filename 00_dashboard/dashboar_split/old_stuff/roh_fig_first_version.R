fractions <- bind_rows(
  
  # compute total fractions
  roh_dat %>% 
    group_by(pop,sample) %>% 
    summarise(per_sample_genome_fraction = sum(length_bp)/genome_length) %>% 
    ungroup() %>% 
    group_by(pop) %>% 
    summarise(mean_per_sample_genome_fraction = mean(per_sample_genome_fraction),
              sd_per_sample_genome_fraction = sd(per_sample_genome_fraction)) %>% 
    mutate(data_set = "total"),
  
  # compute fractions only for ROHs > 2Mbp
  roh_dat %>% 
    filter(length_bp > 1e6) %>% 
    group_by(pop,sample) %>% 
    summarise(per_sample_genome_fraction = sum(length_bp)/genome_length) %>% 
    ungroup() %>% 
    group_by(pop) %>% 
    summarise(mean_per_sample_genome_fraction = mean(per_sample_genome_fraction),
              sd_per_sample_genome_fraction = sd(per_sample_genome_fraction)) %>% 
    mutate(data_set = "at_least_1Mbp"),
  
  # compute fractions only for ROHs > 2Mbp
  roh_dat %>% 
    filter(length_bp > 2e6) %>% 
    group_by(pop,sample) %>% 
    summarise(per_sample_genome_fraction = sum(length_bp)/genome_length) %>% 
    ungroup() %>% 
    group_by(pop) %>% 
    summarise(mean_per_sample_genome_fraction = mean(per_sample_genome_fraction),
              sd_per_sample_genome_fraction = sd(per_sample_genome_fraction)) %>% 
    mutate(data_set = "at_least_2Mbp"),
  
  # compute fractions only for ROHs > 3Mbp
  roh_dat %>% 
    filter(length_bp > 3e6) %>% 
    group_by(pop,sample) %>% 
    summarise(per_sample_genome_fraction = sum(length_bp)/genome_length) %>% 
    ungroup() %>% 
    group_by(pop) %>% 
    summarise(mean_per_sample_genome_fraction = mean(per_sample_genome_fraction),
              sd_per_sample_genome_fraction = sd(per_sample_genome_fraction)) %>% 
    mutate(data_set = "at_least_3Mbp")
  
) %>% 
  mutate(data_set = factor(data_set, levels = c("total","at_least_1Mbp","at_least_2Mbp","at_least_3Mbp"))) %>% 
  mutate(errors_lower = mean_per_sample_genome_fraction - sd_per_sample_genome_fraction,
         errors_upper = mean_per_sample_genome_fraction + sd_per_sample_genome_fraction)


png(here("00_dashboard/for_manuscript/barplot_roh_genomic_fractions.png"), res = 300, units = "px", height = 2000, width = 2000)
ggplot(fractions, aes(x = pop, y = mean_per_sample_genome_fraction, fill = data_set)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = errors_lower, ymax = errors_upper), position =  position_dodge2(width = 0.5, padding = 0.5)) +
  theme_bw(base_size = 12) +
  xlab(NULL) +
  ylab("genome fraction") +
  ggtitle("Mean per-sample genome fractions as ROHs") +
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous(breaks = seq(0,1,0.05), expand = c(0.01,0.01)) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
dev.off()
