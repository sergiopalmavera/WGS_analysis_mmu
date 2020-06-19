# libraries
library(dplyr)
library(stringr)

# Load data
fls <- list.files("../output", pattern = ".pi$")

pi_dat <- lapply(fls, function(fl){
  read.table(file.path("../output", fl), stringsAsFactors = F, header = T)
})

names(pi_dat) <- str_remove(
  fls, "cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.allrecords."
  ) %>% 
  str_remove(".windowed.pi")

# function to prepare data and calculate pi-ratio 
prepare_contrasts <- function(targ_line, ref_line){
#targ_line = "DUK"
#ref_line = "FZTDU"
  
  # merge data
  dd <- full_join(
    pi_dat[[targ_line]],
    pi_dat[[ref_line]],
    by = c("CHROM","BIN_START", "BIN_END"),
    suffix = c(
      paste0("_",targ_line),
      paste0("_",ref_line)
      )
  ) 
  
  # remove entries where reference line (fztdu) is NA
  col_ref <- paste0("PI_",ref_line)
  dd <- dd[ !is.na(dd[,col_ref]) ,]  
  
  # replace NAs in target population with 0 (NA = same DNA in windows)
  col_targ <- paste0("PI_",targ_line)
  dd[ is.na(dd[,col_targ]) , ] <- 0
  
  # choose entries where either of the two lines has >=10 SNPs in window
  idx <- dd[, paste0("N_VARIANTS_",targ_line)] >= 10 | dd[, paste0("N_VARIANTS_",ref_line)] >= 10
  dd <- dd[idx,]
  
  # calculate pi-ratio
  dd$pi_ratio <- dd[,col_targ] / dd[,col_ref]
  
  return(dd)
  }

# prepare pi-ratio tables
pi_ratio_duk_fztdu <- prepare_contrasts("DUK","FZTDU")
pi_ratio_duc_fztdu <- prepare_contrasts("DUC","FZTDU")
pi_ratio_du6_fztdu <- prepare_contrasts("DU6","FZTDU")
pi_ratio_du6p_fztdu <- prepare_contrasts("DU6P","FZTDU")
pi_ratio_duhlb_fztdu <- prepare_contrasts("DUhLB","FZTDU")

# export tables for dashboard
write.csv(pi_ratio_duk_fztdu, "../output/pi_ratio_duk_fztdu.csv")
write.csv(pi_ratio_duc_fztdu, "../output/pi_ratio_duc_fztdu.csv")
write.csv(pi_ratio_du6_fztdu, "../output/pi_ratio_du6_fztdu.csv")
write.csv(pi_ratio_du6p_fztdu, "../output/pi_ratio_du6p_fztdu.csv")
write.csv(pi_ratio_duhlb_fztdu, "../output/pi_ratio_duhlb_fztdu.csv")

#------------------------------------------------------------
# same as in dasboard, but discarded (mostly zeroes!!!)
# i should've expected it, since most sites in selected lines
# are homozygous, while most sites in fztdu are polymorphic
#------------------------------------------------------------

# load data
fls <- list.files("../output", pattern = ".csv$", full.names = T)

pi_ratio_dat <- lapply(fls, function(fl){
  read.csv(fl) %>% dplyr::select(CHROM, BIN_START, BIN_END, pi_ratio)
})

names(pi_ratio_dat) <- basename(fls) %>% str_remove("pi_ratio_") %>% str_remove(".csv")

pi_ratio_dat <- pi_ratio_dat %>% 
  bind_rows(.id = "contrast") %>% 
  mutate(pi_ratio_log2 = log2(pi_ratio))

# summary
summary_win_pi_ratio <- pi_ratio_dat %>% 
  mutate(is_x = (CHROM == "X")) %>% 
  group_by(contrast, is_x) %>% 
  summarise(
    n_win = n(),
    min = min(pi_ratio),
    q1 = quantile(pi_ratio, 0.01),
    q5 = quantile(pi_ratio, 0.05),
    median = median(pi_ratio),
    mean = mean(pi_ratio),
    max = max(pi_ratio)
  ) %>% 
  arrange(is_x)

summary_win_pi_ratio %>%
  mutate(is_x = ifelse(is_x, "X", "Autosome")) %>% 
  dplyr::rename(chr_type = is_x) %>% 
  kable(digits = 2, caption = "Mean Window pi_ratio")


# distribution
p1 <- pi_ratio_dat %>% 
  ggplot(aes(x = pi_ratio)) +
  geom_histogram(bins = 50) +
  facet_wrap(~contrast, ncol = 1, strip.position = "right") +
  theme_bw() +
  xlab("pi-ratio")

p2 <- pi_ratio_dat %>% 
  ggplot(aes(x = NA, y = pi_ratio)) +
  geom_boxplot(width = 0.2) +
  coord_flip() +
  facet_wrap(~contrast, ncol = 1, strip.position = "right") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  ylab("pi-ratio")

gridExtra::grid.arrange(p1,p2, nrow = 1)

rm(p1,p2)


