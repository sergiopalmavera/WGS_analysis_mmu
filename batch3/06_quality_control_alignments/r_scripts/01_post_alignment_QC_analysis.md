Quality Control Alignments
================
Sergio E. Palma-Vera

Load libraries
==============

``` r
library(dplyr)
library(here)
library(stringr)
library(ggplot2)
library(reshape2)
library(tidyr)
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

Load CollectWgsMetrics summary
==============================

``` r
fls <- list.files(here("batch3/06_quality_control_alignments/output"), pattern = "/*_pt1_points.txt")

cvg_summary <- lapply(fls, function(fl){
  
#fl=fls[1]
  # Prepare sample name (taking into account droput-reseqd sample)
  cond <- fl == "I34772-L1_S19_L004.sorted.RG.dedup.bqsr.CollectWgsMetrics.default_pt1_points.txt" | fl == "I34772-L1_S19_L004.sorted.RG.dedup.bqsr.CollectWgsMetrics.Q0.M0_pt1_points.txt"
  if(cond){
    s <- "I34772_dropout_reseqd"
  }else{
    s <- str_split(fl,"-") %>% sapply(function(x) x[1])
  }

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
dim(cvg_summary) #182 31
```

    ## [1] 182  31

Add mouse line group to drop out sample
=======================================

``` r
cvg_summary <- cvg_summary %>% 
  mutate(Linie = ifelse(sample_name == "I34772_dropout_reseqd", "DUK", Linie))
```

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

Summarize (excluding original dropout sample)

``` r
cvg_summary %>% 
  filter(sample_name != "I34772") %>% 
  group_by(mode) %>% 
  summarise(n(), min(MEAN_COVERAGE), mean(MEAN_COVERAGE), median(MEAN_COVERAGE),max(MEAN_COVERAGE),
            sum(MEAN_COVERAGE < 5)) %>% knitr::kable()
```

| mode    |  n()|  min(MEAN\_COVERAGE)|  mean(MEAN\_COVERAGE)|  median(MEAN\_COVERAGE)|  max(MEAN\_COVERAGE)|  sum(MEAN\_COVERAGE &lt; 5)|
|:--------|----:|--------------------:|---------------------:|-----------------------:|--------------------:|---------------------------:|
| Q0.M0   |   90|             2.808120|              7.888959|                8.021949|             13.87711|                          14|
| Q20.M20 |   90|             2.457122|              7.168783|                7.285225|             12.80411|                          19|

``` r
cvg_summary %>% 
  filter(sample_name != "I34772") %>% 
  group_by(Linie,mode) %>% 
  summarise(n(), min(MEAN_COVERAGE), mean(MEAN_COVERAGE), median(MEAN_COVERAGE),max(MEAN_COVERAGE),
            sum(MEAN_COVERAGE < 5)) %>% knitr::kable()
```

| Linie | mode    |  n()|  min(MEAN\_COVERAGE)|  mean(MEAN\_COVERAGE)|  median(MEAN\_COVERAGE)|  max(MEAN\_COVERAGE)|  sum(MEAN\_COVERAGE &lt; 5)|
|:------|:--------|----:|--------------------:|---------------------:|-----------------------:|--------------------:|---------------------------:|
| DU6   | Q0.M0   |   15|             2.808120|              4.491793|                4.613744|             5.704492|                          10|
| DU6   | Q20.M20 |   15|             2.457122|              4.086525|                4.210245|             5.195715|                          13|
| DU6P  | Q0.M0   |   15|             4.004283|              8.758193|                9.103018|            12.336217|                           1|
| DU6P  | Q20.M20 |   15|             3.533915|              7.931205|                8.274174|            11.192941|                           2|
| DUC   | Q0.M0   |   15|             3.915564|              7.422962|                7.215993|            10.871049|                           1|
| DUC   | Q20.M20 |   15|             3.577197|              6.755091|                6.587519|             9.924003|                           2|
| DUHLB | Q0.M0   |   15|             6.590221|              8.273129|                7.808289|            11.819198|                           0|
| DUHLB | Q20.M20 |   15|             5.905461|              7.507391|                7.117772|            10.782619|                           0|
| DUK   | Q0.M0   |   15|             3.549499|              8.041708|                8.182620|            10.690693|                           2|
| DUK   | Q20.M20 |   15|             3.188033|              7.286432|                7.349724|             9.701134|                           2|
| FZTDU | Q0.M0   |   15|             7.490512|             10.345966|                9.922551|            13.877113|                           0|
| FZTDU | Q20.M20 |   15|             6.875381|              9.446055|                9.109785|            12.804107|                           0|

``` r
cvg_summary %>% filter(grepl("I34772", sample_name)) %>% 
  dplyr::select(Linie, sample_name, MEAN_COVERAGE)
```

    ##   Linie           sample_name MEAN_COVERAGE
    ## 1   DUK I34772_dropout_reseqd      3.188033
    ## 2   DUK I34772_dropout_reseqd      3.549499
    ## 3   DUK                I34772      0.183444
    ## 4   DUK                I34772      0.208006

-   In average the mean cvg per sample was between 7-8x (median is similar)

-   However there is a number of samples with less than 5x avg-cvg (14 if M0-Q0, and 19 if M20-Q20)

-   The line with the most number of samples below avg-cvg 5x was DU6 (10 if M0-Q0 and 13 if M20-Q20).

-   The drop-out sample has avg-cvg of ~0.2x, after resequencing mean cvg post-alignment was less than 4 (3,18-3.54).

-   Samples in FZTDU and DUHLB were consistently above 5x

Remove dropout sample
=====================

``` r
cvg_summary <- cvg_summary %>% filter(sample_name != "I34772")
```

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

Remove dropout sample before sequencing

``` r
avg_cvg <- avg_cvg %>% 
  filter(sample_name != "I34772_dropout") 
```

Calculate global cvg (all samples in cvg subset) mean

``` r
global_cvg <- avg_cvg %>% 
  group_by(cvg_class) %>% 
  summarise(mean = mean(avg_cvg))
global_cvg %>% knitr::kable()
```

| cvg\_class      |      mean|
|:----------------|---------:|
| raw             |  9.028242|
| corrected       |  8.734409|
| aligned\_Q0M0   |  7.888959|
| aligned\_Q20M20 |  7.168783|

Organize group levels

``` r
avg_cvg$Linie <- factor(avg_cvg$Linie, levels = c("DUK","DUC","DU6","DU6P","DUHLB","FZTDU"))
```

Visualize

``` r
png(here("batch3/06_quality_control_alignments/figures/avg_cvg_by_class.png"),res = 300, units = "px", height = 2000, width = 3000)

ggplot(data = avg_cvg, aes(x = Linie, y = avg_cvg)) +
    geom_violin() +
    geom_jitter(width = 0.05, alpha = 0.5) +
    geom_boxplot(alpha = 0.2, width = 0.3, outlier.shape = NA) +
    facet_wrap(~cvg_class, nrow=1) +
    geom_hline(data = global_cvg, aes(yintercept = mean), 
               color = "darkgreen", linetype = "dotdash", size = 1) +
    geom_hline(yintercept = 5, color = "red", linetype = "dotted", size = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = seq(0,20,1)) +
    xlab(NULL)

dev.off()
```

-   The mean coverage across all samples ranges from ~9x when data is raw tp 7.1x when reads had been aligned and filter by quality.

-   Read-cleaning (correction) has a slight effect on global avg-cvg (9.02 vs 8.73)

-   Using only high quality base pairs (M20Q20) reduces the aligned global avg-cvg slightly from 7.88x to 7.16x.

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
cvg_summary$Linie <- factor(cvg_summary$Linie, 
                            levels = c("DUK","DUC","DU6","DU6P","DUHLB","FZTDU"))

cvg_summary %>% 
  dplyr::select(mode, Linie, sample_name, cols) %>%
  melt(id.vars = c("mode", "Linie","sample_name")) %>%
  filter(variable == "PCT_5X") %>% 
  group_by(mode, Linie) %>% 
  summarise(n=n(), min=min(value), mean=mean(value), median=median(value), max=max(value)) %>% 
  knitr::kable(digits = 2)
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(cols)` instead of `cols` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

| mode    | Linie |    n|   min|  mean|  median|   max|
|:--------|:------|----:|-----:|-----:|-------:|-----:|
| Q0.M0   | DUK   |   15|  0.32|  0.74|    0.80|  0.89|
| Q0.M0   | DUC   |   15|  0.36|  0.75|    0.78|  0.92|
| Q0.M0   | DU6   |   15|  0.20|  0.45|    0.47|  0.64|
| Q0.M0   | DU6P  |   15|  0.37|  0.79|    0.85|  0.92|
| Q0.M0   | DUHLB |   15|  0.60|  0.76|    0.76|  0.91|
| Q0.M0   | FZTDU |   15|  0.64|  0.82|    0.82|  0.95|
| Q20.M20 | DUK   |   15|  0.28|  0.69|    0.75|  0.85|
| Q20.M20 | DUC   |   15|  0.33|  0.70|    0.73|  0.88|
| Q20.M20 | DU6   |   15|  0.16|  0.40|    0.43|  0.59|
| Q20.M20 | DU6P  |   15|  0.32|  0.74|    0.80|  0.88|
| Q20.M20 | DUHLB |   15|  0.55|  0.72|    0.72|  0.87|
| Q20.M20 | FZTDU |   15|  0.60|  0.78|    0.77|  0.92|

``` r
cvg_summary %>% 
  dplyr::select(mode, Linie, sample_name, cols) %>%
  melt(id.vars = c("mode", "Linie","sample_name")) %>%
  filter(variable == "PCT_5X") %>% 
  group_by(mode) %>% 
  summarise(n=n(), min=min(value), mean=mean(value), median=median(value), max=max(value)) %>% 
  knitr::kable(digits = 2)
```

| mode    |    n|   min|  mean|  median|   max|
|:--------|----:|-----:|-----:|-------:|-----:|
| Q0.M0   |   90|  0.20|  0.72|    0.77|  0.95|
| Q20.M20 |   90|  0.16|  0.67|    0.72|  0.92|

``` r
png(here("batch3/06_quality_control_alignments/figures/pct_territory_min_5x.png"),res = 300, units = "px", height = 2000, width = 3000)

cvg_summary %>% 
  dplyr::select(mode, Linie, sample_name, cols) %>%
  melt(id.vars = c("mode", "Linie","sample_name")) %>%
  filter(variable == "PCT_5X") %>% 
  ggplot(aes(Linie, value*100)) +
    geom_violin() +
    geom_jitter(width = 0.05, alpha = 0.5) +
    geom_boxplot(alpha = 0.2, width = 0.3, outlier.shape = NA) +
    facet_wrap(~mode) +
    ylab("Percentage Genome") +
    xlab(NULL) +
    ggtitle("Percentage Genome Covered at least 5x")

dev.off()
```

-   The average proportion of genome covered at least 5x in any sample was approx 70%.

-   Half of the genomes had more than 70% 5x-cvg.

-   A few genomes had very good cvg with &gt; 90% territory with at least 5x

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
    ## 1       0.967033     85.7

``` r
pct_prop_paired %>% filter(pct == min(pct_prop_paired$pct))
```

    ##   Linie sample_name  pct
    ## 1  DU6P      I34741 85.7

-   96% of samples had a proportion of properly paired reads &gt; 90%

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
  
    cond <- fl == "I34772-L1_S19_L004.sorted.RG.dedup.bqsr.CollectInsertSizeMetrics_pt1.txt"
  if(cond){
    s <- "I34772_dropout_reseqd"
  }else{
    s <- str_split(fl,"-") %>% sapply(function(x) x[1])
  }

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

    ## [1] 91 25

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

-   Most samples had average insert size of ~400bp

-   This is a reasonable interval consistent with the library preparation

Conclusions
===========

-   In average the mean cvg per sample was between 7-8x (median is similar). However there is a number of samples with less than 5x avg-cvg (14 if M0-Q0, and 19 if M20-Q20). The line with the most number of samples below avg-cvg 5x was DU6 (10 if M0-Q0 and 13 if M20-Q20). The drop-out sample has avg-cvg of ~0.2x, after resequencing mean cvg post-alignment was less than 4 (3.18-3.54). Samples in FZTDU and DUHLB were consistently above 5x

-   The mean coverage across all samples ranges from ~9x when data is raw tp 7.1x when reads had been aligned and filter by quality. Correction/cleaning has a slight effect on global avg-cvg (9.02 vs 8.73). Using only high quality base pairs (M20Q20) reduces the aligned global avg-cvg slightly from 7.88x to 7.16x.

-   The average proportion of genome covered at least 5x in any sample was approx 70%. Half of the genomes had more than 70% 5x-cvg. A few genomes had very good cvg with &gt; 90% territory with at least 5x.

-   90% of corrected reads were properly in &gt;95% of the samples. The lowest sample was I34741 (85.7%). Still an acceptable level. The drop-out sample had at least good mapping pct.

-   Most samples had average insert size of ~400bp. This is a reasonable interval consistent with the library preparation
