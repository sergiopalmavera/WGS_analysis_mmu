Quality Control Reads
================
Sergio E. Palma-Vera

Intro
=====

This analysis processes FastQC results on reads before and after quality-trimming and adapter-removal with fastp.

The aim of this analysis is to have an overview of the read quality and coverage (before alignment).

Resequenced-dropout sample (I34772) is included in the data set

Load packages
=============

``` r
library(dplyr)
library(ggplot2)
library(fastqcr)
library(here)
library(reshape2)
library(stringr)
library(pheatmap)
```

Load sample information
=======================

``` r
sample_info <- read.csv(here("sample_info","sample_info_batch3","sample_info.csv"), header = T, stringsAsFactors = F) %>% mutate(name_short = str_remove(name, "-S1")) 

sample_info$Linie[sample_info$Linie == "HLB"] <- "DUHLB"

line_lev <- c("DUK","DUC","DU6","DU6P","DUHLB","FZTDU")

sample_info %>% head()
```

    ##   X      name barcode    plate row col Linie lfd..Nr.  Lab_ID Gen mittl.Verw
    ## 1 1 I34710-S1 S145092 PS000610   A   1 FZTDU        2  FZTDU2  30 0.07277285
    ## 2 2 I34711-S1 S145093 PS000610   B   1 FZTDU        6  FZTDU6  30 0.07221722
    ## 3 3 I34712-S1 S145094 PS000610   C   1 FZTDU        8  FZTDU8  30 0.07146059
    ## 4 4 I34713-S1 S145095 PS000610   D   1 FZTDU       11 FZTDU11  30 0.07003995
    ## 5 5 I34714-S1 S145096 PS000610   E   1 FZTDU       12 FZTDU12  30 0.07219722
    ## 6 6 I34715-S1 S145097 PS000610   F   1 FZTDU       13 FZTDU13  30 0.07187205
    ##   Extraction_date   DNA     Unit X260.280 X260.230 Sample.Type order Plate
    ## 1      04.12.2019  67.5 ng/\xb5l     1.88     2.37         DNA     1 PS611
    ## 2      04.12.2019  89.8 ng/\xb5l     1.86     2.09         DNA     2 PS611
    ## 3      04.12.2019 123.8 ng/\xb5l     1.89     2.37         DNA     3 PS611
    ## 4      04.12.2019  84.8 ng/\xb5l     1.86     2.41         DNA     4 PS611
    ## 5      04.12.2019 100.0 ng/\xb5l     1.87     2.00         DNA     5 PS611
    ## 6      04.12.2019 125.2 ng/\xb5l     1.91     2.43         DNA     6 PS611
    ##   name_short
    ## 1     I34710
    ## 2     I34711
    ## 3     I34712
    ## 4     I34713
    ## 5     I34714
    ## 6     I34715

``` r
sample_info %>% dim()
```

    ## [1] 90 20

Genome mmu length (golden path)
===============================

``` r
genome_length <-    2730871774
```

<http://www.ensembl.org/Mus_musculus/Info/Annotation>

Import fastqc data
==================

Raw reads
---------

``` r
qc_raw_aggr <- qc_aggregate(here("batch3","01_quality_control","output"), progressbar = F)
qc_raw_stats <- qc_stats(qc_raw_aggr)
```

Post-fastp reads
----------------

``` r
qc_fastp_aggr <- qc_aggregate(here("batch3","02_quality_trimming_adapter_removal","output_fastqc"), progressbar = F)
qc_fastp_stats <- qc_stats(qc_fastp_aggr)
```

Total number of files
=====================

``` r
qc_raw_stats %>% nrow()
```

    ## [1] 182

``` r
qc_fastp_stats %>% nrow()
```

    ## [1] 182

There are 182 fastq files because the one of the samples (I34772) was resequenced because it was a drop-out, thus 90+1 \* 2

Count number of files that passed and failed at each module
===========================================================

``` r
summary(qc_raw_aggr)
```

    ## # A tibble: 11 x 7
    ## # Groups:   module [11]
    ##    module    nb_samples nb_fail nb_pass nb_warn failed           warned         
    ##    <chr>          <dbl>   <dbl>   <dbl>   <dbl> <chr>            <chr>          
    ##  1 Adapter …        182       0     182       0 <NA>             <NA>           
    ##  2 Basic St…        182       0     182       0 <NA>             <NA>           
    ##  3 Kmer Con…        182      51     130       1 I34710-L1_S1_L0… I34772-L1_S63_…
    ##  4 Overrepr…        182       8      89      85 I34725-L1_S16_L… I34710-L1_S1_L…
    ##  5 Per base…        182       0     182       0 <NA>             <NA>           
    ##  6 Per base…        182       0     182       0 <NA>             <NA>           
    ##  7 Per base…        182       0     182       0 <NA>             <NA>           
    ##  8 Per sequ…        182     144      33       5 I34711-L1_S2_L0… I34710-L1_S1_L…
    ##  9 Per sequ…        182       0     182       0 <NA>             <NA>           
    ## 10 Sequence…        182       0     182       0 <NA>             <NA>           
    ## 11 Sequence…        182       0     182       0 <NA>             <NA>

Kmer amd per-sequence GC content are FAIL. Kmer content is usually not informative, most people ingore it when it fails. For "per-sequence GC content" there is a second shoulder we have already observed, but nothing major to worry about.

``` r
summary(qc_fastp_aggr)
```

    ## # A tibble: 11 x 7
    ## # Groups:   module [11]
    ##    module    nb_samples nb_fail nb_pass nb_warn failed           warned         
    ##    <chr>          <dbl>   <dbl>   <dbl>   <dbl> <chr>            <chr>          
    ##  1 Adapter …        182       0     182       0 <NA>             <NA>           
    ##  2 Basic St…        182       0     182       0 <NA>             <NA>           
    ##  3 Kmer Con…        182      27     153       2 I34719-L1_S10_L… I34772-L1_S63_…
    ##  4 Overrepr…        182       0     169      13 <NA>             I34725-L1_S16_…
    ##  5 Per base…        182       0     182       0 <NA>             <NA>           
    ##  6 Per base…        182       0     182       0 <NA>             <NA>           
    ##  7 Per base…        182       0     182       0 <NA>             <NA>           
    ##  8 Per sequ…        182     150      29       3 I34711-L1_S2_L0… I34713-L1_S4_L…
    ##  9 Per sequ…        182       0     182       0 <NA>             <NA>           
    ## 10 Sequence…        182       0     182       0 <NA>             <NA>           
    ## 11 Sequence…        182       0       0     182 <NA>             I34710-L1_S1_L…

Kmer content FAIL cases are reduced by half and "per-sequence GC content" FAIL cases increased by 4 cases. But as mentioned earlier, this is an issue we have to live with and it does not mean there is contamination.

Visualize files that passed and failed at each module
=====================================================

Recode status (pass, warn, fail) into integers (2,1,0)
------------------------------------------------------

Raw reads

``` r
qc_raw_aggr$status_recoded[qc_raw_aggr$status == "PASS"] <- 2L
```

    ## Warning: Unknown or uninitialised column: 'status_recoded'.

``` r
qc_raw_aggr$status_recoded[qc_raw_aggr$status == "WARN"] <- 1L
qc_raw_aggr$status_recoded[qc_raw_aggr$status == "FAIL"] <- 0L

qc_raw_aggr_wide <- qc_raw_aggr %>% 
  dplyr::select(module,sample,status_recoded) %>% 
  dcast(module ~ sample, value.var = "status_recoded")

qc_raw_aggr_wide %>% dim()
```

    ## [1]  11 183

``` r
qc_raw_aggr_wide[1:5,1:5]
```

    ##                      module I34710-L1_S1_L001_R1_001 I34710-L1_S1_L001_R2_001
    ## 1           Adapter Content                        2                        2
    ## 2          Basic Statistics                        2                        2
    ## 3              Kmer Content                        0                        2
    ## 4 Overrepresented sequences                        2                        1
    ## 5        Per base N content                        2                        2
    ##   I34711-L1_S2_L001_R1_001 I34711-L1_S2_L001_R2_001
    ## 1                        2                        2
    ## 2                        2                        2
    ## 3                        2                        2
    ## 4                        2                        1
    ## 5                        2                        2

Post fastp reads

``` r
qc_fastp_aggr$status_recoded[qc_fastp_aggr$status == "PASS"] <- 2L
```

    ## Warning: Unknown or uninitialised column: 'status_recoded'.

``` r
qc_fastp_aggr$status_recoded[qc_fastp_aggr$status == "WARN"] <- 1L
qc_fastp_aggr$status_recoded[qc_fastp_aggr$status == "FAIL"] <- 0L

qc_fastp_aggr_wide <- qc_fastp_aggr %>% 
  dplyr::select(module,sample,status_recoded) %>% 
  dcast(module ~ sample, value.var = "status_recoded")

qc_fastp_aggr_wide %>% dim()
```

    ## [1]  11 183

``` r
qc_fastp_aggr_wide[1:5,1:5]
```

    ##                      module I34710-L1_S1_L001_R1_001.corrected
    ## 1           Adapter Content                                  2
    ## 2          Basic Statistics                                  2
    ## 3              Kmer Content                                  2
    ## 4 Overrepresented sequences                                  2
    ## 5        Per base N content                                  2
    ##   I34710-L1_S1_L001_R2_001.corrected I34711-L1_S2_L001_R1_001.corrected
    ## 1                                  2                                  2
    ## 2                                  2                                  2
    ## 3                                  2                                  2
    ## 4                                  2                                  2
    ## 5                                  2                                  2
    ##   I34711-L1_S2_L001_R2_001.corrected
    ## 1                                  2
    ## 2                                  2
    ## 3                                  2
    ## 4                                  2
    ## 5                                  2

Combine raw and post-fastp reads
--------------------------------

``` r
qc_module_status_ma <- left_join(qc_raw_aggr_wide, qc_fastp_aggr_wide, by = "module")

qc_module_status_ma %>% dim()
```

    ## [1]  11 365

Add rownames
------------

``` r
rownames(qc_module_status_ma) <- qc_module_status_ma$module %>% str_replace_all(pattern = " ","_")

qc_module_status_ma$module <- NULL

qc_module_status_ma %>% dim()
```

    ## [1]  11 364

Prepare column meta information
-------------------------------

Get sample name

``` r
samp_nm <- names(qc_module_status_ma) %>% str_split("-") %>% sapply(function(x) x[1])
```

Get read info

``` r
read_info <- names(qc_module_status_ma) %>% str_split("_") %>% sapply(function(x) x[4])
```

Get "stage" info (pre/post fastp)

``` r
qc_stage <- names(qc_module_status_ma) %>% str_split("[.]") %>% sapply(function(x) x[2])

qc_stage[is.na(qc_stage)] <- "raw"

qc_stage[qc_stage == "corrected"] <- "after_fastp"
```

``` r
col_data <- data.frame(read_info = read_info, qc_stage = qc_stage)

rownames(col_data) <- names(qc_module_status_ma)
```

Visualize in a heatmap
----------------------

all 360 files (90+1 samples \* 2 stages \* 2 reads = 360)

``` r
png(here("batch3/03_quality_control_analysis/figures/modules_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma,
         main = "Module status all fastq files\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = col_data, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

raw files (90+1 samples \* 1 stages \* 2 reads = 180)

``` r
rownm_tmp <- rownames(col_data)[col_data$qc_stage == "raw"] # need to adjust coldata
tmp <- filter(col_data, qc_stage == "raw") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_raw_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,-grep("corrected",names(qc_module_status_ma))],
         main = "Module status raw fastq files\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

corrected files (90+1 samples \* 1 stages \* 2 reads = 180)

``` r
rownm_tmp <- rownames(col_data)[col_data$qc_stage == "after_fastp"] # need to adjust coldata
tmp <- filter(col_data, qc_stage == "after_fastp") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_corrected_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,grep("corrected",names(qc_module_status_ma))],
         main = "Module status corrected fastq files\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

all R1 files

``` r
rownm_tmp <- rownames(col_data)[col_data$read_info == "R1"] # need to adjust coldata
tmp <- filter(col_data, read_info == "R1") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_R1_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,grep("R1",names(qc_module_status_ma))],
         main = "Module status R1 fastq files\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

all R2 files

``` r
rownm_tmp <- rownames(col_data)[col_data$read_info == "R2"] # need to adjust coldata
tmp <- filter(col_data, read_info == "R2") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_R2_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,grep("R2",names(qc_module_status_ma))],
         main = "Module status R2 fastq files\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

Prepare qc\_stats
=================

Prepare raw data
----------------

``` r
qc_raw_stats <- qc_raw_stats %>% 
  mutate(tot.seq = as.numeric(tot.seq),
         seq.length = as.numeric(seq.length),
         tot_bp = tot.seq * seq.length,
         read_type = str_split(sample,"_") %>% sapply(function(x) x[4]),
         sample_name = str_split(sample,"_") %>% sapply(function(x) x[1]),
         sample_name = str_remove(sample_name, "-L1")) %>% 
  left_join(dplyr::select(sample_info, name_short, Linie), 
            by = c("sample_name"="name_short")) %>% 
  mutate(Linie = factor(Linie, levels = line_lev))
```

Prepare corrected data
----------------------

Fix columns

``` r
qc_fastp_stats <- qc_fastp_stats %>% 
  mutate(tot.seq = as.numeric(tot.seq),
         seq.length = as.numeric(seq.length),
         tot_bp = tot.seq * seq.length,
         read_type = str_split(sample,"_") %>% sapply(function(x) x[4]),
         sample_name = str_split(sample,"_") %>% sapply(function(x) x[1]),
         sample_name = str_remove(sample_name, "-L1")) %>% 
  left_join(dplyr::select(sample_info, name_short, Linie), 
            by = c("sample_name"="name_short")) %>% 
  mutate(Linie = factor(Linie, levels = line_lev))
```

    ## Warning: NAs durch Umwandlung erzeugt

The columns `seq.length` and `tot_bp` are NAs. The former due to the fact that after fastp reads are of variable sizes - nothing to do about that. However, for the later, I can add the number of bp in file, that have been pre-computed.

Load total base pairs per fastq:

``` r
nbp <- list.files(here("batch3/02_quality_trimming_adapter_removal/output_n_bp"), pattern = ".NBP") %>% 
  lapply(function(x) read.table(file.path(here("batch3/02_quality_trimming_adapter_removal/output_n_bp"),x))) %>% 
  bind_rows()

nbp <- nbp %>% 
  mutate(file = list.files(here("batch3/02_quality_trimming_adapter_removal/output_n_bp"), pattern = ".NBP")) %>% 
  rename(tot_bp = V1) %>% 
  dplyr::select(file, tot_bp) %>% 
  mutate(file = str_remove(file, ".NBP"))
```

Add new total bp info to `qc_fastp_stats`

``` r
qc_fastp_stats <- left_join(qc_fastp_stats, nbp, by = c("sample" = "file"))
qc_fastp_stats <- qc_fastp_stats %>% 
  dplyr::select(-tot_bp.x) %>% 
  rename(tot_bp = tot_bp.y)
head(qc_fastp_stats)
```

    ## # A tibble: 6 x 9
    ##   sample   pct.dup pct.gc tot.seq seq.length read_type sample_name Linie  tot_bp
    ##   <chr>      <dbl>  <dbl>   <dbl>      <dbl> <chr>     <chr>       <fct>   <dbl>
    ## 1 I34710-…    12.0     42  1.31e8         NA R1        I34710      FZTDU 1.96e10
    ## 2 I34710-…    11.6     42  1.31e8         NA R2        I34710      FZTDU 1.96e10
    ## 3 I34711-…    13.7     42  8.80e7         NA R1        I34711      FZTDU 1.32e10
    ## 4 I34711-…    13.0     42  8.80e7         NA R2        I34711      FZTDU 1.32e10
    ## 5 I34712-…    15.7     42  9.19e7         NA R1        I34712      FZTDU 1.37e10
    ## 6 I34712-…    15.0     42  9.19e7         NA R2        I34712      FZTDU 1.37e10

``` r
dim(qc_fastp_stats)
```

    ## [1] 182   9

Label drop out sample, before and after requencing
==================================================

Drop out sample is I34772 (originally as `I34772-L1_S63_L003`, after resequencing `I34772-L1_S19_L004`)

Re-label

``` r
qc_raw_stats %>% filter(sample_name == "I34772")
```

    ## # A tibble: 4 x 9
    ##   sample    pct.dup pct.gc tot.seq seq.length tot_bp read_type sample_name Linie
    ##   <chr>       <dbl>  <dbl>   <dbl>      <dbl>  <dbl> <chr>     <chr>       <fct>
    ## 1 I34772-L…   14.5      43  3.77e7        151 5.69e9 R1        I34772      DUK  
    ## 2 I34772-L…   13.2      43  3.77e7        151 5.69e9 R2        I34772      DUK  
    ## 3 I34772-L…    7.85     43  2.18e6        151 3.29e8 R1        I34772      DUK  
    ## 4 I34772-L…    7.13     43  2.18e6        151 3.29e8 R2        I34772      DUK

``` r
qc_raw_stats <- qc_raw_stats %>% 
    mutate(sample_name = ifelse(sample == "I34772-L1_S19_L004_R1_001","I34772_dropout_reseq", sample_name),
           sample_name = ifelse(sample == "I34772-L1_S19_L004_R2_001", "I34772_dropout_reseq", sample_name),
           sample_name = ifelse(sample == "I34772-L1_S63_L003_R1_001", "I34772_dropout", sample_name),
           sample_name = ifelse(sample == "I34772-L1_S63_L003_R2_001", "I34772_dropout", sample_name))

qc_raw_stats %>% filter(grepl("I34772", sample_name))
```

    ## # A tibble: 4 x 9
    ##   sample   pct.dup pct.gc  tot.seq seq.length tot_bp read_type sample_name Linie
    ##   <chr>      <dbl>  <dbl>    <dbl>      <dbl>  <dbl> <chr>     <chr>       <fct>
    ## 1 I34772-…   14.5      43 37688140        151 5.69e9 R1        I34772_dro… DUK  
    ## 2 I34772-…   13.2      43 37688140        151 5.69e9 R2        I34772_dro… DUK  
    ## 3 I34772-…    7.85     43  2181992        151 3.29e8 R1        I34772_dro… DUK  
    ## 4 I34772-…    7.13     43  2181992        151 3.29e8 R2        I34772_dro… DUK

``` r
qc_fastp_stats %>% filter(sample_name == "I34772")
```

    ## # A tibble: 4 x 9
    ##   sample    pct.dup pct.gc tot.seq seq.length read_type sample_name Linie tot_bp
    ##   <chr>       <dbl>  <dbl>   <dbl>      <dbl> <chr>     <chr>       <fct>  <dbl>
    ## 1 I34772-L…   14.2      43  3.68e7         NA R1        I34772      DUK   5.51e9
    ## 2 I34772-L…   12.8      43  3.68e7         NA R2        I34772      DUK   5.51e9
    ## 3 I34772-L…    7.45     43  2.13e6         NA R1        I34772      DUK   3.18e8
    ## 4 I34772-L…    6.66     43  2.13e6         NA R2        I34772      DUK   3.18e8

``` r
qc_fastp_stats <- qc_fastp_stats %>% 
  mutate(sample_name = ifelse(sample == "I34772-L1_S19_L004_R1_001.corrected", "I34772_dropout_reseq", sample_name),
         sample_name = ifelse(sample == "I34772-L1_S19_L004_R2_001.corrected", "I34772_dropout_reseq", sample_name),
         sample_name = ifelse(sample == "I34772-L1_S63_L003_R1_001.corrected", "I34772_dropout", sample_name),
         sample_name = ifelse(sample == "I34772-L1_S63_L003_R2_001.corrected", "I34772_dropout", sample_name))

qc_fastp_stats %>% filter(grepl("I34772", sample_name))
```

    ## # A tibble: 4 x 9
    ##   sample    pct.dup pct.gc tot.seq seq.length read_type sample_name Linie tot_bp
    ##   <chr>       <dbl>  <dbl>   <dbl>      <dbl> <chr>     <chr>       <fct>  <dbl>
    ## 1 I34772-L…   14.2      43  3.68e7         NA R1        I34772_dro… DUK   5.51e9
    ## 2 I34772-L…   12.8      43  3.68e7         NA R2        I34772_dro… DUK   5.51e9
    ## 3 I34772-L…    7.45     43  2.13e6         NA R1        I34772_dro… DUK   3.18e8
    ## 4 I34772-L…    6.66     43  2.13e6         NA R2        I34772_dro… DUK   3.18e8

Evaluate number of pairs and cvg per sample from raw fastq files
================================================================

Same number of reads in R1 and R2 files?

``` r
qc_raw_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot.seq") %>% 
  apply(1, function(x) x[3]==x[4]) %>% all()
```

    ## [1] TRUE

Yes

Visualize number of read pairs per sample

``` r
png("./batch3/03_quality_control_analysis/figures/n_read_pairs_raw_fastq.png", 
    height = 2500, width = 4000, units = "px", res = 300)
qc_raw_stats %>% 
  filter(read_type == "R1") %>% 
  ggplot(aes(x = sample_name, y = tot.seq/1e6, fill = Linie)) + 
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("N Read Pairs (in mill)") +
    ggtitle("Number of read pairs in million in raw fastq files")
dev.off()
```

Visualize number of reads per sample

``` r
png("./batch3/03_quality_control_analysis/figures/n_reads_raw_fastq.png", 
    height = 2500, width = 4000, units = "px", res = 300)
qc_raw_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot.seq") %>%
  mutate(tot_reads_samp = R1 + R2) %>% 
  ggplot(aes(x = sample_name, y = tot_reads_samp, fill = Linie)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("total read per sample (in mill)") +
    ggtitle("Number of reads in million in raw fastq files")
dev.off()
```

-   One of the samples (I34772, DUK) has much less number of pair reads than the rest of the samples.

Visualize coverage

``` r
png("./batch3/03_quality_control_analysis/figures/avg_cvg_raw_fastq.png", 
    height = 2500, width = 4000, units = "px", res = 300)
qc_raw_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  ggplot(aes(x = sample_name, y = avg_cvg, fill = Linie)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Average Genome Coverage per sample") +
    ggtitle("Average Genome Coverage per sample from raw fastq files") +
    geom_hline(yintercept = 5, color = "red", linetype = "dashed")
dev.off()    
```

Number of samples below avg cvg 5x

``` r
qc_raw_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  group_by(Linie) %>% 
  summarise(n_below_5x = sum(avg_cvg < 5))
```

    ## # A tibble: 6 x 2
    ##   Linie n_below_5x
    ##   <fct>      <int>
    ## 1 DUK            2
    ## 2 DUC            1
    ## 3 DU6            4
    ## 4 DU6P           1
    ## 5 DUHLB          0
    ## 6 FZTDU          0

There are 8 samples below the targeted avg-cvg

Check coverage of those samples

``` r
qc_raw_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  filter(avg_cvg < 5)
```

    ##            sample_name Linie         R1         R2      tot_bp   avg_cvg
    ## 1               I34725   DU6 5490718776 5490718776 10981437552 4.0212205
    ## 2               I34727   DU6 6189158706 6189158706 12378317412 4.5327348
    ## 3               I34728   DU6 6723928528 6723928528 13447857056 4.9243825
    ## 4               I34733   DU6 6204795662 6204795662 12409591324 4.5441867
    ## 5               I34753  DU6P 6665741735 6665741735 13331483470 4.8817684
    ## 6               I34764   DUC 6224749255 6224749255 12449498510 4.5588001
    ## 7       I34772_dropout   DUK  329480792  329480792   658961584 0.2413008
    ## 8 I34772_dropout_reseq   DUK 5690909140 5690909140 11381818280 4.1678333

Out of the 8 samples, 7 samples are above 4x, while one is only 0.24x

The drop-out appears here twice as below 5x, after resequencing, the sample remains &lt;5x.

Evaluate number of pairs and cvg per sample from corrected fastq files
======================================================================

Same number of reads in R1 and R2 files?

``` r
qc_fastp_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot.seq") %>% 
  apply(1, function(x) x[3]==x[4]) %>% all()
```

    ## [1] TRUE

Yes

Visualize number of read pairs per sample

``` r
png("./batch3/03_quality_control_analysis/figures/n_read_pairs_corrected_fastq.png", 
    height = 2500, width = 4000, units = "px", res = 300)
qc_fastp_stats %>% 
  filter(read_type == "R1") %>% 
  ggplot(aes(x = sample_name, y = tot.seq/1e6, fill = Linie)) + 
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("N Read Pairs (in mill)") +
    ggtitle("Number of read pairs in million in corrected fastq files")
dev.off()
```

Visualize coverage

``` r
png("./batch3/03_quality_control_analysis/figures/avg_cvg_corrected_fastq.png", 
    height = 2500, width = 4000, units = "px", res = 300)
qc_fastp_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  ggplot(aes(x = sample_name, y = avg_cvg, fill = Linie)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Average Genome Coverage per sample") +
    ggtitle("Average Genome Coverage per sample from corrected fastq files") +
    geom_hline(yintercept = 5, color = "red", linetype = "dashed")
dev.off()    
```

Number of samples below avg cvg 5x

``` r
qc_fastp_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  group_by(Linie) %>% 
  summarise(n_below_5x = sum(avg_cvg < 5))
```

    ## # A tibble: 6 x 2
    ##   Linie n_below_5x
    ##   <fct>      <int>
    ## 1 DUK            2
    ## 2 DUC            1
    ## 3 DU6            5
    ## 4 DU6P           1
    ## 5 DUHLB          0
    ## 6 FZTDU          0

-   There are now 9 samples below the targeted avg-cvg

Check coverage of those samples

``` r
qc_fastp_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  filter(avg_cvg < 5)
```

    ##            sample_name Linie         R1         R2      tot_bp   avg_cvg
    ## 1               I34725   DU6 4883094508 4880278420  9763372928 3.5751854
    ## 2               I34727   DU6 5790293841 5789699950 11579993791 4.2404019
    ## 3               I34728   DU6 6217699942 6216400929 12434100871 4.5531617
    ## 4               I34732   DU6 6557408375 6556425506 13113833881 4.8020687
    ## 5               I34733   DU6 5881036486 5880272312 11761308798 4.3067964
    ## 6               I34753  DU6P 6438166208 6437956573 12876122781 4.7150228
    ## 7               I34764   DUC 6090560999 6090838719 12181399718 4.4606268
    ## 8       I34772_dropout   DUK  318197885  318184757   636382642 0.2330328
    ## 9 I34772_dropout_reseq   DUK 5506818339 5506549391 11013367730 4.0329128

-   Out of the 9 samples, 7 samples are above 4x, one sample is above 3x and one is only 0.23x

-   The drop-out appears here twice as below 5x, after resequencing, the sample remains &lt;5x.

-   The new sample below 5x is I34732

Summarize GC content
====================

``` r
qc_raw_stats %>% summarise(mean(pct.gc))
```

    ## # A tibble: 1 x 1
    ##   `mean(pct.gc)`
    ##            <dbl>
    ## 1           42.3

``` r
qc_fastp_stats %>% summarise(mean(pct.gc))
```

    ## # A tibble: 1 x 1
    ##   `mean(pct.gc)`
    ##            <dbl>
    ## 1           42.1

-   The mean GC content does not change much after fastp-correction.

-   This level of GC content is consistent with the expectation for mouse (42%, <https://bionumbers.hms.harvard.edu/bionumber.aspx>?&id=102409&ver=9)

Summarize pct duplication
=========================

``` r
qc_raw_stats %>% summarise(mean(pct.dup))
```

    ## # A tibble: 1 x 1
    ##   `mean(pct.dup)`
    ##             <dbl>
    ## 1            12.4

``` r
qc_fastp_stats %>% summarise(mean(pct.dup))
```

    ## # A tibble: 1 x 1
    ##   `mean(pct.dup)`
    ##             <dbl>
    ## 1            12.0

-   The duplication pct is ~12%. So probably there will be 12% drop in coverage after alignment.

Conclusions
===========

-   Raw fastq files and corrected fastq files each corresponds to 180 files (90 samples, R1 and R2 reads). All modules were visualized on a heatmap colored according to the module-status (PASS, WARN, FAIL).

-   Overall these modules had the PASS status before and after correction with fasp: "adapter content", "basic statisitcs", "per base n content", "per base sequence content", "per base sequence quality", "per sequence quality scores", "sequence duplication levels"

-   The module "sequence length distribution" changed from PASS to WARN in all samples after fastp-correction. But this is nothing to worry about; it is a consequence of the quality trimming and adapter removal, because reads become shorter.

-   The module "kmer content" was partially fixed after correction.

-   The "overrepresented sequences" WARN/FAIL was partially solved by fastp (only a few WARN remeained).

-   "Per sequence GC content" fails before and after correction, because there is a peak and a "shoulder". This happens in most of the samples. The peak aligns with the expected distribution's peak. I had a similar problem before, and it is not something to worry about (<https://www.biostars.org/p/341611/>). We also discussed this in one of the status meetings and we agreed that it required no further attention. The worst case scenario is that it is due to contamination, but those sequences would no align to the reference. Reads have a mean GC content (as seen in the fastp reports) close to the expected GC content for mouse (42%, <https://bionumbers.hms.harvard.edu/bionumber.aspx>?&id=102409&ver=9). Since this situation is very similar to what was observed before, it requires no further attention.

-   Other than the GC-distribution issue, which requires no further measures, the data seems good enough to proceed.

-   In terms of raw avg-cvg, there were 7 samples below 5x, out of which 6 samples were above 4x and one sample (DUK) only had 0.24x

-   Corrected (fastq files after fastp) avg-cvg showed 8 samples below 5x, 6 were above 4x, one was above 3x and one had only 0.23x. The new sample below 5x is sample I34732 coresponding to a DU6 sample.

-   The mean GC content does not change much after fastp-correction (42.3 vs 42.1). This level of GC content is consistent with the expectation for mouse (42%, <https://bionumbers.hms.harvard.edu/bionumber.aspx>?&id=102409&ver=9)

-   The duplication pct is ~12%. So probably there will be ~12% drop in coverage after alignment.

-   After resequencing the dropout sample (original avg cvg 0.23x) coverage increased to 4x (original dropout not included, as it could be a suboptimal sample).

Export cvg information
======================

``` r
qc_raw_stats_avg_cvg <- qc_raw_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  dplyr::select(-R1,-R2,-tot_bp)
```

``` r
qc_fastp_stats_avg_cvg <- qc_fastp_stats %>% 
  dcast(sample_name+Linie ~ read_type, value.var = "tot_bp") %>% 
  mutate(tot_bp = R1 + R2, avg_cvg = tot_bp/genome_length) %>% 
  dplyr::select(-R1,-R2,-tot_bp)
```

``` r
saveRDS(qc_raw_stats_avg_cvg, here("batch3/03_quality_control_analysis/r_objects/qc_raw_stats_avg_cvg.rds"))
saveRDS(qc_fastp_stats_avg_cvg, here("batch3/03_quality_control_analysis/r_objects/qc_fastp_stats_avg_cvg.rds"))
```
