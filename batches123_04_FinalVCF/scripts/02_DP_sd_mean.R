library(dplyr)
library(vroom)
library(here)

dat <- vroom(here("batches123_04_FinalVCF/output/DP.tmp"), col_names = F)

# last col is not needed
res_mean <- dat[,1:150] %>% as.matrix() %>% as.numeric() %>% mean(na.rm=T)
res_mean

res_sd <- dat[,1:150] %>% as.matrix() %>% as.numeric() %>% sd(na.rm=T)
res_sd

# Max DP according to https://www.biostars.org/p/265782/
max_DP <- res_mean + 4*res_sd
max_DP




