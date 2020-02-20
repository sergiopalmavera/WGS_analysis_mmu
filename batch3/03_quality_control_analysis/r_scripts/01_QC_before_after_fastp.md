Quality Control Reads
================
Sergio E. Palma-Vera

Intro
=====

This analysis processes FastQC results on reads before and after quality-trimming and adapter-removal with fastp.

The aim of this analysis is to have an overview of the read quality and coverage (before alignment).

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
sample_info <- read.csv(here("sample_info","sample_info_batch3","sample_info.csv"), header = T, stringsAsFactors = F)
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

``` r
sample_info %>% dim()
```

    ## [1] 90 19

Import fastqc data
==================

Raw reads

``` r
qc_raw_aggr <- qc_aggregate(here("batch3","01_quality_control","output"), progressbar = F)
qc_raw_stats <- qc_stats(qc_raw_aggr)
```

Post-fastp reads

``` r
qc_fastp_aggr <- qc_aggregate(here("batch3","02_quality_trimming_adapter_removal","output_fastqc"), progressbar = F)
qc_fastp_stats <- qc_stats(qc_fastp_aggr)
```

Count number of files that passed and failed at each module
===========================================================

``` r
summary(qc_raw_aggr)
```

    ## # A tibble: 11 x 7
    ## # Groups:   module [11]
    ##    module    nb_samples nb_fail nb_pass nb_warn failed           warned         
    ##    <chr>          <dbl>   <dbl>   <dbl>   <dbl> <chr>            <chr>          
    ##  1 Adapter …        180       0     180       0 <NA>             <NA>           
    ##  2 Basic St…        180       0     180       0 <NA>             <NA>           
    ##  3 Kmer Con…        180      51     128       1 I34710-L1_S1_L0… I34772-L1_S63_…
    ##  4 Overrepr…        180       8      88      84 I34725-L1_S16_L… I34710-L1_S1_L…
    ##  5 Per base…        180       0     180       0 <NA>             <NA>           
    ##  6 Per base…        180       0     180       0 <NA>             <NA>           
    ##  7 Per base…        180       0     180       0 <NA>             <NA>           
    ##  8 Per sequ…        180     142      33       5 I34711-L1_S2_L0… I34710-L1_S1_L…
    ##  9 Per sequ…        180       0     180       0 <NA>             <NA>           
    ## 10 Sequence…        180       0     180       0 <NA>             <NA>           
    ## 11 Sequence…        180       0     180       0 <NA>             <NA>

``` r
summary(qc_fastp_aggr)
```

    ## # A tibble: 11 x 7
    ## # Groups:   module [11]
    ##    module    nb_samples nb_fail nb_pass nb_warn failed           warned         
    ##    <chr>          <dbl>   <dbl>   <dbl>   <dbl> <chr>            <chr>          
    ##  1 Adapter …        180       0     180       0 <NA>             <NA>           
    ##  2 Basic St…        180       0     180       0 <NA>             <NA>           
    ##  3 Kmer Con…        180      27     151       2 I34719-L1_S10_L… I34772-L1_S63_…
    ##  4 Overrepr…        180       0     167      13 <NA>             I34725-L1_S16_…
    ##  5 Per base…        180       0     180       0 <NA>             <NA>           
    ##  6 Per base…        180       0     180       0 <NA>             <NA>           
    ##  7 Per base…        180       0     180       0 <NA>             <NA>           
    ##  8 Per sequ…        180     148      29       3 I34711-L1_S2_L0… I34713-L1_S4_L…
    ##  9 Per sequ…        180       0     180       0 <NA>             <NA>           
    ## 10 Sequence…        180       0     180       0 <NA>             <NA>           
    ## 11 Sequence…        180       0       0     180 <NA>             I34710-L1_S1_L…

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

    ## [1]  11 181

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

    ## [1]  11 181

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

    ## [1]  11 361

Add rownames
------------

``` r
rownames(qc_module_status_ma) <- qc_module_status_ma$module %>% str_replace_all(pattern = " ","_")

qc_module_status_ma$module <- NULL

qc_module_status_ma %>% dim()
```

    ## [1]  11 360

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

Visualize in a heatmap all 360 files (90 samples \* 2 stages \* 2 reads = 360)
------------------------------------------------------------------------------

``` r
png(here("batch3/03_quality_control_analysis/figures/modules_360_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma,
         main = "Module status all fastq files (n=360)\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = col_data, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

    ## png 
    ##   3

Visualize in a heatmap raw files (90 samples \* 1 stages \* 2 reads = 180)
--------------------------------------------------------------------------

``` r
rownm_tmp <- rownames(col_data)[col_data$qc_stage == "raw"] # need to adjust coldata
tmp <- filter(col_data, qc_stage == "raw") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_180_raw_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,-grep("corrected",names(qc_module_status_ma))],
         main = "Module status raw fastq files (n=180)\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

Visualize in a heatmap corrected files (90 samples \* 1 stages \* 2 reads = 180)
--------------------------------------------------------------------------------

``` r
rownm_tmp <- rownames(col_data)[col_data$qc_stage == "after_fastp"] # need to adjust coldata
tmp <- filter(col_data, qc_stage == "after_fastp") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_180_corrected_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,grep("corrected",names(qc_module_status_ma))],
         main = "Module status corrected fastq files (n=180)\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

Visualize all R1 files
----------------------

``` r
rownm_tmp <- rownames(col_data)[col_data$read_info == "R1"] # need to adjust coldata
tmp <- filter(col_data, read_info == "R1") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_180_R1_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,grep("R1",names(qc_module_status_ma))],
         main = "Module status R1 fastq files (n=180)\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

Visualize all R2 files
----------------------

``` r
rownm_tmp <- rownames(col_data)[col_data$read_info == "R2"] # need to adjust coldata
tmp <- filter(col_data, read_info == "R2") 
rownames(tmp) <- rownm_tmp

png(here("batch3/03_quality_control_analysis/figures/modules_180_R2_files.png"), 
    width = 2000, height = 1500, units = "px", res = 300)
pheatmap(qc_module_status_ma[,grep("R2",names(qc_module_status_ma))],
         main = "Module status R2 fastq files (n=180)\n(PASS=2 / WARN=1 / FAIL=0)",
         annotation_col = tmp, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F)
dev.off()
```

Conclusions
===========

-   Raw fastq files and corrected fastq files each corresponds to 180 files. All modules were visualized on a heatmap colored according to the module-status (PASS, WARN, FAIL).

-   Overall these modules had the PASS status before and after correction with fasp: "adapter content", "basic statisitcs", "per base n content", "per base sequence content", "per base sequence quality", "per sequence quality scores", "sequence duplication levels"

-   The module "sequence length distribution" changed from PASS to WARN in all samples after fastp-correction. But this is nothing to worry about; it is a consequence of the quality trimming and adapter removal, because reads become shorter.

-   The module "kmer content" was partially fixed after correction.

-   The "overrepresented sequences" WARN/FAIL was partially solved by fastp (only a few WARN remeained).

-   "Per sequence GC content" fails before and after correction, because there is a peak and a "shoulder". This happens in most of the samples. The peak aligns with the expected distribution's peak. I had a similar problem before, and it is not something to worry about (<https://www.biostars.org/p/341611/>). We also discussed this in one of the status meetings and we agreed that it required no further attention. The worst case scenario is that it is due to contamination, but those sequences would no align to the reference. Reads have a mean GC content (as seen in the fastp reports) close to the expected GC content for mouse (42%, <https://bionumbers.hms.harvard.edu/bionumber.aspx>?&id=102409&ver=9). Since this situation is very similar to what was observed before, it requires no further attention.

-   Other than the GC-distribution issue, which requires no further measures, the data seems good enough to proceed.
