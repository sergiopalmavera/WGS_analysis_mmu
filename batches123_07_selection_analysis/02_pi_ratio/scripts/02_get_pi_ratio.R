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

# export tables
write.csv(pi_ratio_duk_fztdu, "../output/pi_ratio_duk_fztdu.csv")
write.csv(pi_ratio_duc_fztdu, "../output/pi_ratio_duc_fztdu.csv")
write.csv(pi_ratio_du6_fztdu, "../output/pi_ratio_du6_fztdu.csv")
write.csv(pi_ratio_du6p_fztdu, "../output/pi_ratio_du6p_fztdu.csv")
write.csv(pi_ratio_duhlb_fztdu, "../output/pi_ratio_duhlb_fztdu.csv")





