library(vroom)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)

# name of data file
#tab <- "tst.filtered_vcf.table"
tab <- "../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.table"

# sample info batch1-2 (60 samples, ~20x)
sample_info1 <- read.csv("../../sample_info/sample_info_batch1_batch2.csv", stringsAsFactors = F)

# sample info last batch (90 samples, ~5x avg-cvg)
sample_info2 <- read.csv("../../sample_info/sample_info_batch3/sample_info.csv", stringsAsFactors = F)

# Combine both data sets
ss1 <- sample_info1 %>% 
  dplyr::select(Linie, sample_id) %>% 
  mutate(target_cvg = "30x")

ss2 <- sample_info2 %>% dplyr::select(Linie, name) %>% 
  mutate(name = str_remove(name,"-S1")) %>% 
  dplyr::rename(sample_id = name) %>% 
  mutate(target_cvg = "5x",
         Linie = ifelse(Linie == "HLB", "DUhLB", Linie)) 

sample_info <- bind_rows(ss1,ss2) %>% 
  mutate(Linie = factor(Linie, c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU"))) %>% 
  arrange(Linie, sample_id) %>% 
  mutate(Linie = as.character(Linie))

dim(sample_info)

# Load genotype table
gt_tab <- vroom(tab)

# Prepare col names
names(gt_tab) <- str_remove(names(gt_tab), ".GT") %>% str_remove("-L1")

# Match order of table and sample info
idx1 <- match(sample_info$sample_id, names(gt_tab))
tmp1 <- gt_tab[,idx1]
gt_tab <- gt_tab %>% 
  dplyr::select(CHROM,POS) %>% 
  bind_cols(tmp1) 

# columns in the right order?
identical(sample_info$sample_id,names(tmp1))

# make a function to extract number of missing samples per pop at each SNP
get_n_miss <- function(l){
  i_col <- names(gt_tab) %in% sample_info$sample_id[sample_info$Linie == l]
  gt_tab[,i_col] %>% apply(1,function(row_x) sum(row_x == "./."))
}

# Put count missings per populations counts into one data frame
miss_cts <- data.frame(
  DUK = get_n_miss("DUK"),
  DUC = get_n_miss("DUC"),
  DU6 = get_n_miss("DU6"),
  DU6P = get_n_miss("DU6P"),
  DUhLB = get_n_miss("DUhLB"),
  FZTDU = get_n_miss("FZTDU")
) 

miss_cts2 <- miss_cts %>%
  # melt for plotting
  melt(variable.name = "Line", value.name = "miss_counts") %>%
  mutate(Line = factor(Line, levels = c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU") ))

png("../figures_tables/miss_cts.png", width = 3000, height = 2000, units = "px", res = 300)
miss_cts2 %>%
  # Visualize
  ggplot(aes(x = miss_counts, after_stat(count), color = Line)) +
    geom_density() +
    #facet_wrap(~Line, ncol = 1, strip.position = "right") +
    theme_grey(base_size = 15) +
    scale_x_continuous(breaks = seq(0,30,1)) +
    xlab("missing counts") +
    ylab("SNP counts") +
    ggtitle("Distribution of Number of Samples Missing per SNP")
dev.off()    


png("../figures_tables/miss_cts_facet.png", width = 3000, height = 3000, units = "px", res = 300)
miss_cts2 %>%
  # Visualize
  ggplot(aes(x = miss_counts)) +
	geom_histogram(binwidth=1) +
	facet_wrap(~Line, ncol = 1, strip.position = "right") +
	theme_grey(base_size = 15) +
	scale_x_continuous(breaks = seq(0,30,1)) +
	xlab("missing counts") +
	ylab("SNP counts") +
	ggtitle("Distribution of Number of Samples Missing per SNP")
dev.off() 
