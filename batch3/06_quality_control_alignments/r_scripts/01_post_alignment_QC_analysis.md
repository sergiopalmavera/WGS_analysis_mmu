Quality Control Alignments
================
Sergio E. Palma-Vera

Load libraries
==============

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(here)
```

    ## here() starts at /projekte/I2-SOS-FERT/GitHub/WGS_analysis_mmu

``` r
library(stringr)
library(ggplot2)
library(reshape2)
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:reshape2':
    ## 
    ##     smiths

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

Load CollectWgsMetrics summary
==============================

``` r
fls <- list.files(here("batch3/06_quality_control_alignments/output"), pattern = "/*_pt1_points.txt")

cvg_summary <- lapply(fls, function(fl){
  
#fl=fls[1]

  # Prepare sample name
  s <- str_split(fl,"-") %>% sapply(function(x) x[1])
  
  # Load sample histogram information
  d <- read.table(file.path(here("batch3/06_quality_control_alignments/output"),fl), header = T, stringsAsFactors = F)
  
  # Extract mode of metrics (with or without bp and mapping quality filters)
  m <- str_split(fl,"_") %>% sapply(function(x) x[3]) %>% 
    str_remove("sorted.RG.dedup.bqsr.CollectWgsMetrics.") %>% 
    str_remove("L00[1-9].")
  
  # Add sample and mode information (default = mapping and bp quality of min 20)
  d <- mutate(d, sample = s, mode = m, mode = str_replace(mode, "default","Q20.M20"))
  
  # add mouse line information
  d <- dplyr::select(sample_info, Linie, name_short) %>% 
    right_join(d, by = c("name_short"="sample")) %>% 
    dplyr::rename(sample_name = name_short)
  
}) %>% 
  # bind rows
  bind_rows()

head(cvg_summary)
```

    ##   Linie sample_name GENOME_TERRITORY MEAN_COVERAGE SD_COVERAGE MEDIAN_COVERAGE
    ## 1 FZTDU      I34710       2652783500     12.131920    6.546794              12
    ## 2 FZTDU      I34710       2652783500     13.192371    6.969047              13
    ## 3 FZTDU      I34711       2652783500      7.859407    5.537658               7
    ## 4 FZTDU      I34711       2652783500      8.623020    6.054945               8
    ## 5 FZTDU      I34712       2652783500      8.010361    6.289709               7
    ## 6 FZTDU      I34712       2652783500      8.813493    6.859339               8
    ##   MAD_COVERAGE PCT_EXC_MAPQ PCT_EXC_DUPE PCT_EXC_UNPAIRED PCT_EXC_BASEQ
    ## 1            4     0.081131     0.049572         0.000365      0.023780
    ## 2            4     0.000000     0.071091         0.000459      0.000019
    ## 3            4     0.105121     0.050443         0.000362      0.023950
    ## 4            4     0.000000     0.085824         0.000464      0.000008
    ## 5            4     0.122014     0.053757         0.000324      0.022503
    ## 6            4     0.000000     0.102177         0.000417      0.000008
    ##   PCT_EXC_OVERLAP PCT_EXC_CAPPED PCT_EXC_TOTAL   PCT_1X   PCT_5X  PCT_10X
    ## 1        0.013445       0.002717      0.171009 0.934378 0.889265 0.678519
    ## 2        0.016463       0.009880      0.097912 0.962141 0.927665 0.731288
    ## 3        0.014115       0.004221      0.198213 0.922374 0.709245 0.356582
    ## 4        0.018031       0.015012      0.119340 0.952303 0.754760 0.402089
    ## 5        0.011757       0.004805      0.215160 0.908281 0.662648 0.366522
    ## 6        0.015469       0.017066      0.135138 0.938751 0.706914 0.409334
    ##    PCT_15X  PCT_20X  PCT_25X  PCT_30X  PCT_40X  PCT_50X  PCT_60X  PCT_70X
    ## 1 0.347626 0.105744 0.019548 0.003481 0.000978 0.000592 0.000403 0.000297
    ## 2 0.399164 0.135958 0.031462 0.008250 0.002546 0.001356 0.000874 0.000635
    ## 3 0.106849 0.017437 0.002530 0.000984 0.000509 0.000321 0.000230 0.000177
    ## 4 0.134518 0.027949 0.006363 0.002834 0.001199 0.000693 0.000478 0.000362
    ## 5 0.151248 0.041051 0.007658 0.001602 0.000529 0.000330 0.000229 0.000169
    ## 6 0.181407 0.056239 0.013832 0.004338 0.001450 0.000773 0.000504 0.000367
    ##    PCT_80X  PCT_90X PCT_100X HET_SNP_SENSITIVITY HET_SNP_Q    mode
    ## 1 0.000232 0.000188 0.000154            0.896523        10 Q20.M20
    ## 2 0.000489 0.000392 0.000323            0.928478        11   Q0.M0
    ## 3 0.000143 0.000120 0.000103            0.776194         7 Q20.M20
    ## 4 0.000287 0.000235 0.000198            0.808792         7   Q0.M0
    ## 5 0.000129 0.000105 0.000088            0.736411         6 Q20.M20
    ## 6 0.000283 0.000229 0.000191            0.770044         6   Q0.M0

``` r
dim(cvg_summary) #180 31
```

    ## [1] 180  31

What's the mean coverage per sample?
====================================

Mean coverage per sample as boxplot/violinplot

``` r
png(here("batch3/06_quality_control_alignments/figures/boxplot_mean_coverage_per_sample.png"),res = 300, units = "px", height = 2000, width = 2000)

cvg_summary %>% 
  ggplot(data = ., aes(x = mode, y = MEAN_COVERAGE)) +
    geom_violin() +
    geom_jitter(width = 0.05, alpha = 0.5, aes(color = Linie )) +
    geom_boxplot(alpha = 0.2, width = 0.05, outlier.shape = NA) +
    xlab(NULL) +
    geom_hline(yintercept = 5, color = "red", linetype = "dotted")

dev.off()
```

Mean coverage per sample as barplot

``` r
png(here("batch3/06_quality_control_alignments/figures/barplot_mean_coverage_per_sample.png"),res = 300, units = "px", height = 2000, width = 3500)

cvg_summary %>% 
  ggplot(data = ., aes(x = sample_name, y = MEAN_COVERAGE, fill = Linie)) +
    xlab(NULL) +
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 5, color = "red", linetype = "dotted") +
    facet_wrap(~mode, nrow = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust  = 0.5))

dev.off()
```

Summarize

``` r
cvg_summary %>% 
  group_by(mode) %>% 
  summarise(n(), min(MEAN_COVERAGE), mean(MEAN_COVERAGE), median(MEAN_COVERAGE),max(MEAN_COVERAGE),
            sum(MEAN_COVERAGE < 5))
```

    ## # A tibble: 2 x 7
    ##   mode  `n()` `min(MEAN_COVER… `mean(MEAN_COVE… `median(MEAN_CO…
    ##   <chr> <int>            <dbl>            <dbl>            <dbl>
    ## 1 Q0.M0    90            0.208             7.85             8.02
    ## 2 Q20.…    90            0.183             7.14             7.29
    ## # … with 2 more variables: `max(MEAN_COVERAGE)` <dbl>, `sum(MEAN_COVERAGE <
    ## #   5)` <int>

-   In average the mean cvg per sample was between 7-8x (median is similar)

-   However there is a number of samples with less than 5x avg-cvg (14 if M0-Q0, and 19 if M20-Q20)

-   The line with the most number of samples below avg-cvg 5x was DU6 (10 if M0-Q0 and 13 if M20-Q20).

-   The drop-out sample has avg-cvg of ~0.2x

-   Samples in FZTDU and DUHLB were consistently above 5x

Compare mean cvg raw reads, corrected reads and aligned reads
=============================================================

Load coverage information before alignment (raw fastq and fastp-corrected fastqs) and combine avg-cvg's

``` r
avg_cvg <- list(
  aligned_Q0M0 = cvg_summary %>% 
    filter(mode == "Q0.M0") %>% 
    dplyr::select(sample_name,Linie, MEAN_COVERAGE) %>% 
    dplyr::rename(avg_cvg = MEAN_COVERAGE),
  aligned_Q20M20 = cvg_summary %>% 
    filter(mode == "Q20.M20") %>% 
    dplyr::select(sample_name,Linie, MEAN_COVERAGE) %>% 
    dplyr::rename(avg_cvg = MEAN_COVERAGE),  
  corrected = readRDS(here("batch3/03_quality_control_analysis/r_objects/qc_fastp_stats_avg_cvg.rds")),
  raw = readRDS(here("batch3/03_quality_control_analysis/r_objects/qc_raw_stats_avg_cvg.rds"))
  ) %>% 
  bind_rows(.id = "cvg_class") %>% 
  mutate(cvg_class = factor(cvg_class, levels = c("raw","corrected","aligned_Q0M0","aligned_Q20M20")))
```

    ## Warning in bind_rows_(x, .id): binding character and factor vector, coercing
    ## into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector, coercing
    ## into character vector

``` r
head(avg_cvg)  
```

    ##      cvg_class sample_name Linie   avg_cvg
    ## 1 aligned_Q0M0      I34710 FZTDU 13.192371
    ## 2 aligned_Q0M0      I34711 FZTDU  8.623020
    ## 3 aligned_Q0M0      I34712 FZTDU  8.813493
    ## 4 aligned_Q0M0      I34713 FZTDU 13.572355
    ## 5 aligned_Q0M0      I34714 FZTDU  9.922551
    ## 6 aligned_Q0M0      I34715 FZTDU 11.394611

Calculate global (all samples in cvg subset) mean

``` r
global_cvg <- avg_cvg %>% 
  group_by(cvg_class) %>% 
  summarise(mean = mean(avg_cvg))
global_cvg
```

    ## # A tibble: 4 x 2
    ##   cvg_class       mean
    ##   <fct>          <dbl>
    ## 1 raw             8.98
    ## 2 corrected       8.69
    ## 3 aligned_Q0M0    7.85
    ## 4 aligned_Q20M20  7.14

Visualize

``` r
png(here("batch3/06_quality_control_alignments/figures/avg_cvg_by_class.png"),res = 300, units = "px", height = 2000, width = 2000)

ggplot(data = avg_cvg, aes(x = Linie, y = avg_cvg)) +
    geom_violin() +
    geom_jitter(width = 0.05, alpha = 0.5) +
    geom_boxplot(alpha = 0.2, width = 0.3, outlier.shape = NA) +
    facet_wrap(~cvg_class, nrow=1) +
    geom_hline(data = global_cvg, aes(yintercept = mean), color = "red", linetype = "dotted") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(NULL)

dev.off()
```

    ## png 
    ##   2

-   The mean coverage across all samples ranges from ~8.7x when data is raw tp 7.14x when reads had been aligned and filter by quality.

-   Correction has little effect on global avg-cvg (8.98 vs 8.69)

-   Using only high quality base pairs (M20Q20) reduces the aligned global avg-cvg slightly from 7.85x to 7.14x.

Inspect the fraction of bases that attained at least X sequence coverage
========================================================================

``` r
cols <- names(cvg_summary)[grep("^PCT_[[:digit:]]",names(cvg_summary))] # get cols with cvg
```

``` r
bks <- str_remove(cols, "PCT_") %>% str_remove("X")

png(here("batch3/06_quality_control_alignments/figures/min_cvg_genome_pct.png"),res = 300, units = "px", height = 2000, width = 2000)

cvg_summary %>% 
  dplyr::select(mode, Linie, sample_name, cols) %>% 
  melt(id.vars = c("mode", "Linie","sample_name")) %>% 
  mutate(value = as.numeric(value),
         cvg = str_remove(as.character(variable), "PCT_"),
         cvg = str_remove(cvg, "X"),
         cvg = as.numeric(cvg)) %>% 
  filter(cvg <= 20) %>% 
  ggplot(aes(x = cvg, y = value, color = sample_name, group = sample_name)) +
    geom_line() +
    facet_grid(Linie ~ mode) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 5, color = "red", linetype = "dotdash") +
    scale_x_continuous(breaks = 1:20)

dev.off()
```

Faction of genome covered by at least 5x

``` r
cvg_summary %>% 
  dplyr::select(mode, Linie, sample_name, cols) %>%
  melt(id.vars = c("mode", "Linie","sample_name")) %>%
  filter(variable == "PCT_5X") %>% 
  group_by(mode) %>% 
  summarise(n(), min(value), mean(value), median(value), max(value))
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(cols)` instead of `cols` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

    ## # A tibble: 2 x 6
    ##   mode    `n()` `min(value)` `mean(value)` `median(value)` `max(value)`
    ##   <chr>   <int>        <dbl>         <dbl>           <dbl>        <dbl>
    ## 1 Q0.M0      90     0.000163         0.715           0.769        0.951
    ## 2 Q20.M20    90     0.000078         0.668           0.718        0.915

-   The average proportion of genome covered at least 5x in each sample was approx 70%.

-   Half of the genomes had more than 70% 5x-cvg.

-   A few genomes had very good cvg with &gt; 90% territory with at least 5x

Prepare histogram data
======================

Define file names

``` r
fls <- list.files(here("batch3/06_quality_control_alignments/output"), pattern = "CollectWgsMetrics") %>% 
  .[grep("_pt2",.)]
```

Files were extracted from second part of picard's CollectWgsMetrics' output.

Import histogram data and add mouse-line information and cummulative mass

``` r
hist_dat <- lapply(fls, function(fl){
  
  # Prepare sample name
  s <- str_split(fl,"-") %>% sapply(function(x) x[1])
  
  # Load sample histogram information
  d <- read.table(file.path(here("batch3/06_quality_control_alignments/output"),fl), header = T)
  
  # Extract mode of metrics (with or without bp and mapping quality filters)
  m <- str_split(fl,"_") %>% sapply(function(x) x[3]) %>% 
    str_remove("sorted.RG.dedup.bqsr.CollectWgsMetrics.") %>% 
    str_remove("L00[1-9].")
  
  # Add sample and mode information (default = mapping and bp quality of min 20)
  d <- mutate(d, sample = s, mode = m, mode = str_replace(mode, "default","Q20.M20"))
  
  # add mouse line information
  d <- dplyr::select(sample_info, Linie, name_short) %>% 
    right_join(d, by = c("name_short"="sample")) %>% 
    dplyr::rename(sample_name = name_short)
  
  # Add probability and cummulative mass
  d <- d %>% mutate(prob = high_quality_coverage_count/sum(high_quality_coverage_count),
                    cum_mass = cumsum(prob))
  
  # change count variable name
  d <- d %>% rename(cvg_count = high_quality_coverage_count)
}) %>% 
  # bind rows
  bind_rows()

dim(hist_dat) #45180     6
```

    ## [1] 45180     7

``` r
head(hist_dat)
```

    ##   Linie sample_name coverage cvg_count    mode        prob   cum_mass
    ## 1 FZTDU      I34710        0 174079653 Q20.M20 0.065621508 0.06562151
    ## 2 FZTDU      I34710        1  16951031 Q20.M20 0.006389904 0.07201141
    ## 3 FZTDU      I34710        2  21272548 Q20.M20 0.008018954 0.08003037
    ## 4 FZTDU      I34710        3  32532477 Q20.M20 0.012263525 0.09229389
    ## 5 FZTDU      I34710        4  48920890 Q20.M20 0.018441343 0.11073523
    ## 6 FZTDU      I34710        5  68984881 Q20.M20 0.026004716 0.13673995

Visualize cummulative mass
==========================

``` r
png(here("batch3/06_quality_control_alignments/figures/cum_sum_cvg_prob.png"),res = 300, units = "px", height = 2000, width = 2000)

hist_dat %>% 
  filter(coverage <= 50) %>% 
  ggplot(data = ., aes(x = coverage, y = cum_mass, color = sample_name)) +
    geom_line() +
    facet_grid(Linie~mode) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 5, color = "red", linetype = "dotdash")

dev.off()
```

-   All samples had some fraction of their genomes with a coverage less than 5x (the target coverage).

-   For some samples this fraction was as high as 100% (the dropout sample)

-   The mouse line DU6 had all samples with 50% of their genomes covered by less than 5x.

-   For the rest of the lines, almost all genomes had &lt;50% below 5x.

Whats the percentage of the genome that is covered 0x?
======================================================

Visualize distribution of pct of genome with zero cvg

``` r
png(here("batch3/06_quality_control_alignments/figures/violin_zero_cvg.png"),res = 300, units = "px", height = 2000, width = 2000)

hist_dat %>% 
  group_by(mode, sample_name) %>% 
  filter(coverage == 0) %>% 
  mutate(pct_zero_cvg = prob*100) %>% 
  ggplot(aes(x=mode, y = pct_zero_cvg)) +
    geom_violin() +
    geom_jitter(width = 0.05, alpha = 0.5) +
    ggtitle("Percentage of genome with 0x covg")

dev.off()  
```

Summarize pct of genome regions with zero cvg

``` r
hist_dat %>% 
  group_by(mode) %>% 
  filter(coverage == 0) %>% 
  mutate(pct_zero_cvg = prob*100) %>% 
  summarise(min(pct_zero_cvg), median(pct_zero_cvg), mean(pct_zero_cvg), max(pct_zero_cvg))
```

    ## # A tibble: 2 x 5
    ##   mode   `min(pct_zero_cv… `median(pct_zero_… `mean(pct_zero_c… `max(pct_zero_c…
    ##   <chr>              <dbl>              <dbl>             <dbl>            <dbl>
    ## 1 Q0.M0               3.72               4.58              6.38             82.3
    ## 2 Q20.M…              6.47               7.70              9.52             84.0

-   The average pct of genome with zero coverage is between 6.38 -9.52 (median between 4.58 - 7.70). So in average, the proportion of the genomes not covered (0x) is less than 10%.

-   There are a few samples with more than 10% of their genome uncovered.

-   There is one sample (the "drop-out" under re-sequencing) that has &gt; 80% of its genome without coverage.

Whats the percentage of the genome that is covered 5x (the agreed target coverage) or more?
===========================================================================================

Visualize distribution of pct of genome with 5x cvg

``` r
png(here("batch3/06_quality_control_alignments/figures/violin_above_5x_cvg.png"),res = 300, units = "px", height = 2000, width = 2000)

hist_dat %>% 
  group_by(mode, sample_name) %>% 
  filter(coverage >= 5) %>% 
  summarise(pct_above_5x_cvg = sum(prob)*100) %>% 
  ggplot(aes(x=mode, y = pct_above_5x_cvg)) +
    geom_violin() +
    geom_jitter(width = 0.05, alpha = 0.5) +
    ggtitle("Percentage of genome covered by at least 5x")
  
dev.off()  
```

Summarize pct of genome regions with zero cvg

``` r
hist_dat %>% 
  group_by(mode, sample_name) %>% 
  filter(coverage >= 5) %>% 
  summarise(pct_above_5x_cvg = sum(prob)*100)  %>% 
  summarise(n(),min(pct_above_5x_cvg), median(pct_above_5x_cvg), mean(pct_above_5x_cvg), max(pct_above_5x_cvg))
```

    ## # A tibble: 2 x 6
    ##   mode  `n()` `min(pct_above_… `median(pct_abo… `mean(pct_above…
    ##   <chr> <int>            <dbl>            <dbl>            <dbl>
    ## 1 Q0.M0    90          0.0163              76.9             71.5
    ## 2 Q20.…    90          0.00783             71.8             66.8
    ## # … with 1 more variable: `max(pct_above_5x_cvg)` <dbl>

-   In average, the proportion of samples covered by at least 5x is between 71.5 and 66.8%.

-   In 50% of the samples (median) the proportion of genome covered by at least 5x is above 70%.

-   The dropput sample seems to have 0% of genome covered by at least 5x.

-   In average, approx. 25% of the genome is covered less than 5x.

How many reads were lost after mapping?
=======================================

Prepare data

``` r
pct_prop_paired <- read.table(here("batch3/06_quality_control_alignments/output/flagstat_pct_properly_paired_reads.tab"))

names(pct_prop_paired) <- c("sample_name","pct")

pct_prop_paired <- pct_prop_paired %>% left_join(sample_info, by = c("sample_name"="name_short"))
```

    ## Warning: Column `sample_name`/`name_short` joining factor and character vector,
    ## coercing into character vector

``` r
pct_prop_paired <- pct_prop_paired %>% dplyr::select(Linie,sample_name, pct)

head(pct_prop_paired)
```

    ##   Linie sample_name   pct
    ## 1 FZTDU      I34710 96.75
    ## 2 FZTDU      I34711 95.94
    ## 3 FZTDU      I34712 94.83
    ## 4 FZTDU      I34713 97.75
    ## 5 FZTDU      I34714 96.71
    ## 6 FZTDU      I34715 96.28

Visualize

``` r
png(here("batch3/06_quality_control_alignments/figures/pct_properly_paired_reads.png"),
    res = 300, units = "px", height = 1000, width = 3500)

ggplot(pct_prop_paired, aes(x=sample_name, y = pct, fill = Linie)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 100, linetype = "dotted", color = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(breaks = seq(0,100,5)) +
  ggtitle("Percentage of properly paired reads") +
  xlab(NULL)
dev.off()
```

``` r
pct_prop_paired %>% summarise(mean(pct > 90), min(pct))
```

    ##   mean(pct > 90) min(pct)
    ## 1      0.9666667     85.7

``` r
pct_prop_paired %>% filter(pct == min(pct_prop_paired$pct))
```

    ##   Linie sample_name  pct
    ## 1  DU6P      I34741 85.7

-   90% of corrected reads were properly in &gt;95% of the samples

-   The lowest sample was I34741 (85.7%). Still an acceptable level.

-   The drop-out sample had at least good mapping pct.

What size are the insert (interval between mapped pairs)
========================================================

Load data

``` r
fls <- list.files(here("batch3/06_quality_control_alignments/output"), pattern = "CollectInsertSizeMetrics_pt1.txt")

insert_size_summary <- lapply(fls, function(fl){
  
#fl=fls[1]
  # Prepare sample name
  s <- str_split(fl,"-") %>% sapply(function(x) x[1])
  
  # Load sample histogram information
  d <- read.delim(file.path(here("batch3/06_quality_control_alignments/output"),fl), 
                  header = T, stringsAsFactors = F, comment.char = "#")
  
  # Add sample and mode information (default = mapping and bp quality of min 20)
  d <- mutate(d, sample = s)
  
  # add mouse line information
  d <- dplyr::select(sample_info, Linie, name_short) %>% 
    right_join(d, by = c("name_short"="sample")) %>% 
    dplyr::rename(sample_name = name_short)
  
}) %>% 
  # bind rows
  bind_rows()

head(insert_size_summary)
```

    ##   Linie sample_name MEDIAN_INSERT_SIZE MODE_INSERT_SIZE
    ## 1 FZTDU      I34710                387              365
    ## 2 FZTDU      I34711                384              355
    ## 3 FZTDU      I34712                388              366
    ## 4 FZTDU      I34713                398              379
    ## 5 FZTDU      I34714                394              376
    ## 6 FZTDU      I34715                403              377
    ##   MEDIAN_ABSOLUTE_DEVIATION MIN_INSERT_SIZE MAX_INSERT_SIZE MEAN_INSERT_SIZE
    ## 1                        61               2       189829258         399.0558
    ## 2                        61               2       188005223         396.7593
    ## 3                        59               2       191115076         396.2769
    ## 4                        60               2       184827350         407.0728
    ## 5                        59               2       187344453         403.9806
    ## 6                        61               2       188611108         412.0448
    ##   STANDARD_DEVIATION READ_PAIRS PAIR_ORIENTATION WIDTH_OF_10_PERCENT
    ## 1          105.98335  117982190               FR                  25
    ## 2          108.43712   77226715               FR                  23
    ## 3           99.02429   77869177               FR                  23
    ## 4           99.10979  121573271               FR                  23
    ## 5          102.14334   88715912               FR                  23
    ## 6          102.13256  101329704               FR                  25
    ##   WIDTH_OF_20_PERCENT WIDTH_OF_30_PERCENT WIDTH_OF_40_PERCENT
    ## 1                  47                  71                  97
    ## 2                  47                  71                  97
    ## 3                  45                  67                  93
    ## 4                  47                  69                  95
    ## 5                  45                  69                  93
    ## 6                  47                  71                  97
    ##   WIDTH_OF_50_PERCENT WIDTH_OF_60_PERCENT WIDTH_OF_70_PERCENT
    ## 1                 123                 153                 189
    ## 2                 123                 155                 191
    ## 3                 119                 147                 183
    ## 4                 121                 151                 185
    ## 5                 119                 149                 185
    ## 6                 123                 155                 191
    ##   WIDTH_OF_80_PERCENT WIDTH_OF_90_PERCENT WIDTH_OF_95_PERCENT
    ## 1                 239                 325                 433
    ## 2                 243                 337                 455
    ## 3                 229                 313                 405
    ## 4                 233                 313                 399
    ## 5                 233                 319                 417
    ## 6                 239                 321                 415
    ##   WIDTH_OF_99_PERCENT SAMPLE LIBRARY READ_GROUP
    ## 1                 827     NA      NA         NA
    ## 2                 833     NA      NA         NA
    ## 3                 695     NA      NA         NA
    ## 4                 669     NA      NA         NA
    ## 5                 751     NA      NA         NA
    ## 6                 695     NA      NA         NA

``` r
dim(insert_size_summary) #90 25
```

    ## [1] 90 25

Visualize

``` r
png(here("batch3/06_quality_control_alignments/figures/insert_size_mean.png"),
    res = 300, units = "px", height = 1000, width = 3500)

ggplot(insert_size_summary, aes(x=sample_name, y = MEAN_INSERT_SIZE, fill = Linie)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = MEAN_INSERT_SIZE-STANDARD_DEVIATION, 
                    ymax = MEAN_INSERT_SIZE+STANDARD_DEVIATION), width = 0.2, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) +
  ggtitle("Mean Insert Size per sample (+/- sd)")

dev.off()
```

    ## png 
    ##   2

-   Most samples had average insert size of ~400bp

-   This is a reasonable interval consistent with the library preparation

Conclusions
===========

-   In average the mean cvg per sample was between 7-8x (median is similar). However there is a number of samples with less than 5x avg-cvg (14 if M0-Q0, and 19 if M20-Q20). The line with the most number of samples below avg-cvg 5x was DU6 (10 if M0-Q0 and 13 if M20-Q20). The drop-out sample has avg-cvg of ~0.2x. Samples in FZTDU and DUHLB were consistently above 5x

-   The mean coverage across all samples ranges from ~8.7x when data is raw tp 7.14x when reads had been aligned and filter by quality. Correction has little effect on global avg-cvg (8.98 vs 8.69). Using only high quality base pairs (M20Q20) reduces the aligned global avg-cvg slightly from 7.85x to 7.14x.

-   The average proportion of genome covered at least 5x in each sample was approx 70%. Half of the genomes had more than 70% 5x-cvg. A few genomes had very good cvg with &gt; 90% territory with at least 5x. All samples had some fraction of their genomes with a coverage less than 5x (the target coverage). For some samples this fraction was as high as 100% (the dropout sample)

-   The mouse line DU6 had all samples with 50% of their genomes covered by less than 5x. For the rest of the lines, almost all genomes had &lt;50% below 5x. The average pct of genome with zero coverage is between 6.38 -9.52 (median between 4.58 - 7.70). So in average, the proportion of the genomes not covered (0x) is less than 10%. There are a few samples with more than 10% of their genome uncovered (0x). There is one sample (the "drop-out" under re-sequencing) that has &gt; 80% of its genome without coverage. In average, the proportion of samples covered by at least 5x is between 71.5 and 66.8%. In 50% of the samples (median) the proportion of genome covered by at least 5x is above 70%. The dropput sample seems to have 0% of genome covered by at least 5x. In average, approx. 25% of the genome is covered less than 5x.

-   90% of corrected reads were properly in &gt;95% of the samples. The lowest sample was I34741 (85.7%). Still an acceptable level. The drop-out sample had at least good mapping pct.

-   Most samples had average insert size of ~400bp. This is a reasonable interval consistent with the library preparation
