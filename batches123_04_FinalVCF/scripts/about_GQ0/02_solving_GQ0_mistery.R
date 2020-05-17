
library(here)
library(vroom)
library(dplyr)
library(tidyr)

dat <- vroom(here("batches123_04_FinalVCF/scripts/tmp/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered_H07788.table"))

#dat <- dat %>% head()

# Remove nonvariants
dat <- dat %>% 
  mutate(GT = `H07788-L1.GT`) %>% 
  separate(col = `H07788-L1.GT`, sep = "/", into = c("allele_1","allele_2")) %>% 
  filter(allele_1 != REF | allele_2 != REF)

nrow(dat) #2752822

head(dat)

dat$GT %>% unique()


dat_excl_missing <- dat %>% 
  filter(GT != "./.")

nrow(dat_excl_missing) #2706051, same as picard's TOTAL_SNPS column


dat_excl_missing$`H07788-L1.GQ` %>% summary() # there are no samples below 20!...as expected


dat %>% filter(GT == "./." & `H07788-L1.GQ` == 0) %>% nrow() #3464 as reported by picard
dat %>% filter(`H07788-L1.GQ` == 0) %>% nrow() #3464 as reported by picard


# Conclusion
## picard counts only called&non-hom-ref as total SNPs
## These SNPs are allo minGQ20
## The SNPs below 20 (incl 0), come from nocalls (./.)...just ignore this column.



