sample information preparation
================
Sergio E. Palma-Vera
18 2 2020

Load pkgs
=========

``` r
library(dplyr)
library(here)
```

Intro
=====

This script puts together internal (fbn-Henry) information and ikmb information to make a file of sample information.

Load ikmb sample info
=====================

``` r
info_ikmb <- read.csv(here("sample_info","sample_info_batch3","sample_info_ikmb.csv"), header = T, stringsAsFactors = F)

info_ikmb %>% head()
```

    ##        name barcode    plate row col
    ## 1 I34710-S1 S145092 PS000610   A   1
    ## 2 I34711-S1 S145093 PS000610   B   1
    ## 3 I34712-S1 S145094 PS000610   C   1
    ## 4 I34713-S1 S145095 PS000610   D   1
    ## 5 I34714-S1 S145096 PS000610   E   1
    ## 6 I34715-S1 S145097 PS000610   F   1

Load internal info
==================

``` r
info_fbn <- read.csv2(here("sample_info","sample_info_batch3","Maus_DNA für Sequenzierung20182019.csv"), header = T, stringsAsFactors = F)

info_fbn %>% head()
```

    ##   Linie lfd..Nr.  Lab_ID Gen mittl.Verw Extraction_date   DNA     Unit X260.280
    ## 1 FZTDU        2  FZTDU2  30 0.07277285      04.12.2019  67.5 ng/\xb5l     1.88
    ## 2 FZTDU        6  FZTDU6  30 0.07221722      04.12.2019  89.8 ng/\xb5l     1.86
    ## 3 FZTDU        8  FZTDU8  30 0.07146059      04.12.2019 123.8 ng/\xb5l     1.89
    ## 4 FZTDU       11 FZTDU11  30 0.07003995      04.12.2019  84.8 ng/\xb5l     1.86
    ## 5 FZTDU       12 FZTDU12  30 0.07219722      04.12.2019 100.0 ng/\xb5l     1.87
    ## 6 FZTDU       13 FZTDU13  30 0.07187205      04.12.2019 125.2 ng/\xb5l     1.91
    ##   X260.230 Sample.Type order Plate Position1 Position2
    ## 1     2.37         DNA     1 PS611         A         1
    ## 2     2.09         DNA     2 PS611         B         1
    ## 3     2.37         DNA     3 PS611         C         1
    ## 4     2.41         DNA     4 PS611         D         1
    ## 5     2.00         DNA     5 PS611         E         1
    ## 6     2.43         DNA     6 PS611         F         1

Merge
=====

``` r
sample_info <- left_join(info_ikmb, info_fbn, by = c("row"="Position1", "col"="Position2"))

sample_info %>% head()
```

    ##        name barcode    plate row col Linie lfd..Nr.  Lab_ID Gen mittl.Verw
    ## 1 I34710-S1 S145092 PS000610   A   1 FZTDU        2  FZTDU2  30 0.07277285
    ## 2 I34711-S1 S145093 PS000610   B   1 FZTDU        6  FZTDU6  30 0.07221722
    ## 3 I34712-S1 S145094 PS000610   C   1 FZTDU        8  FZTDU8  30 0.07146059
    ## 4 I34713-S1 S145095 PS000610   D   1 FZTDU       11 FZTDU11  30 0.07003995
    ## 5 I34714-S1 S145096 PS000610   E   1 FZTDU       12 FZTDU12  30 0.07219722
    ## 6 I34715-S1 S145097 PS000610   F   1 FZTDU       13 FZTDU13  30 0.07187205
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

    ## [1] 90 18

Export
======

``` r
#write.csv(sample_info, here("sample_info","sample_info.csv"))
```
