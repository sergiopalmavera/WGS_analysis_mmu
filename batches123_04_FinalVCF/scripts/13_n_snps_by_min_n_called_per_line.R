
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-suse-linux-gnu (64-bit)

R ist freie Software und kommt OHNE JEGLICHE GARANTIE.
Sie sind eingeladen, es unter bestimmten Bedingungen weiter zu verbreiten.
Tippen Sie 'license()' or 'licence()' für Details dazu.

R ist ein Gemeinschaftsprojekt mit vielen Beitragenden.
Tippen Sie 'contributors()' für mehr Information und 'citation()',
um zu erfahren, wie R oder R packages in Publikationen zitiert werden können.

Tippen Sie 'demo()' für einige Demos, 'help()' für on-line Hilfe, oder
'help.start()' für eine HTML Browserschnittstelle zur Hilfe.
Tippen Sie 'q()', um R zu verlassen.

[Vorher gesicherter Workspace wiederhergestellt]

> library(vroom)
> library(dplyr)

Attache Paket: 'dplyr'

The following objects are masked from 'package:stats':

    filter, lag

The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union

> library(ggplot2)
> library(stringr)
> library(reshape2)
> 
> # name of data file
> #tab <- "tst.filtered_vcf.table"
> tab <- "../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.table"
> 
> # sample info batch1-2 (60 samples, ~20x)
> sample_info1 <- read.csv("../../sample_info/sample_info_batch1_batch2.csv", stringsAsFactors = F)
> 
> # sample info last batch (90 samples, ~5x avg-cvg)
> sample_info2 <- read.csv("../../sample_info/sample_info_batch3/sample_info.csv", stringsAsFactors = F)
> 
> # Combine both data sets
> ss1 <- sample_info1 %>% 
+   dplyr::select(Linie, sample_id) %>% 
+   mutate(target_cvg = "30x")
> 
> ss2 <- sample_info2 %>% dplyr::select(Linie, name) %>% 
+   mutate(name = str_remove(name,"-S1")) %>% 
+   dplyr::rename(sample_id = name) %>% 
+   mutate(target_cvg = "5x",
+          Linie = ifelse(Linie == "HLB", "DUhLB", Linie)) 
> 
> sample_info <- bind_rows(ss1,ss2) %>% 
+   mutate(Linie = factor(Linie, c("DUK","DUC","DU6","DU6P","DUhLB","FZTDU"))) %>% 
+   arrange(Linie, sample_id) %>% 
+   mutate(Linie = as.character(Linie))
> 
> dim(sample_info)
[1] 150   3
> 
> # Load genotype table
> gt_tab <- vroom(tab)
Rows: 9,115,052
Columns: 152
Delimiter: "\t"
chr [150]: H07738-L1.GT, H07739.GT, H07740-L1.GT, H07741-L1.GT, H07742-L1.GT, H07743-L1.GT, ...
dbl [  2]: CHROM, POS

Use `spec()` to retrieve the guessed column specification
Pass a specification to the `col_types` argument to quiet this message
> 
> # Prepare col names
> names(gt_tab) <- str_remove(names(gt_tab), ".GT") %>% str_remove("-L1")
> 
> # Match order of table and sample info
> idx1 <- match(sample_info$sample_id, names(gt_tab))
> tmp1 <- gt_tab[,idx1]
> gt_tab <- gt_tab %>% 
+   dplyr::select(CHROM,POS) %>% 
+   bind_cols(tmp1) 
> 
> # columns in the right order?
> identical(sample_info$sample_id,names(tmp1))
[1] TRUE
> 
> # make a function to extract number of missing samples per pop at each SNP
> get_n_not_miss <- function(l){
+   i_col <- names(gt_tab) %in% sample_info$sample_id[sample_info$Linie == l]
+   gt_tab[,i_col] %>% apply(1,function(row_x) sum(row_x != "./."))
+ }
> 
> # Put count missings per populations counts into one data frame
> not_miss_cts <- data.frame(
+   DUK = get_n_not_miss("DUK"),
+   DUC = get_n_not_miss("DUC"),
+   DU6 = get_n_not_miss("DU6"),
+   DU6P = get_n_not_miss("DU6P"),
+   DUhLB = get_n_not_miss("DUhLB"),
+   FZTDU = get_n_not_miss("FZTDU")
+ ) 
> 
> # define min number of called-samples in pop 
> minNs <- 0:25
> 
> # Count number of SNPs per population with at least minN calls
> n_snps <- data.frame(
+   min_called_samples = minNs,
+   DUK = vector(mode = "numeric", length = length(minNs)),
+   DUC = vector(mode = "numeric", length = length(minNs)),
+   DU6 = vector(mode = "numeric", length = length(minNs)),
+   DU6P = vector(mode = "numeric", length = length(minNs)),
+   DUhLB = vector(mode = "numeric", length = length(minNs)),
+   FZTDU = vector(mode = "numeric", length = length(minNs))
+ )
> 
> for(i in 1:length(minNs)){
+   
+   minN <- minNs[i]
+   
+   res <- sapply(not_miss_cts, function(col_i) sum(col_i >= minN) ) 
+   
+   n_snps[i,2:ncol(n_snps)] <- res
+ }
> 
> # convert to fraction
> fraction_snps <- lapply(n_snps[2:ncol(n_snps)], function(x) x/nrow(n_snps)) %>% 
+   bind_cols() %>% 
+   mutate(min_called_samples = minNs) %>% 
+   dplyr::select(min_called_samples, everything())
> 
> 
> write.csv(n_snps, "../figures_tables/n_snps_by_min_n_called_per_line.csv")
> 
> write.csv(fraction_snps, "../figures_tables/fraction_snps_by_min_n_called_per_line.csv")
> 
> proc.time()
       User      System verstrichen 
   2407.410      36.451    2362.358 
