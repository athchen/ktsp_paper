# ktsp_paper

This repository contains all of the data and the code to reproduce the work presented in the paper **A Top Scoring Pairs Classifier for Recent HIV Infections**. 

The session info, including all packages used in the analyses and their respective versions can be found in the [Session Info](#session-info) section of this README.md. All code used to generate the results of the manuscript are mapped in the [section below](#mapping-of-manuscript-content-to-repository-files). A detailed description of the contents of each file in the directory is provided in the [Directory Structure](#directory-structure) section. 

# Directory structure

### **`data_raw`**

Contains `.rda` files with normalized read counts and relevant sample and peptide annotation. 

* **`hiv_data.rda`**: data for the training and testing set.
* **`validation_data.rda`**: data for the validation set.

Within each `.rda` file are give objects:

* `rc_ptid_yrs`: matrix of log10 normalized read counts for each samples (rows are samples, columns are peptides).
* `pt_anno`: annotation matrix for all study participants.
* `sample_anno`: sample specific metadata.
* `pep_anno`: description of the peptides.
* `hiv_ind`: row indices of `pep_anno` indicating which peptides span the HIV genome.

### **`data_processed`**

Generated output used in later parts of the analyses.

* **`kpep_list.rds`**: sorted data frame of the top 10 peptide pairs, along with their optimal cutoffs. Generated in `R/run-ktsp.R`.
* **`slopes.rda`**: regression results obtained using code in `R/run-slopes.R`. 

### **`figures`**

Contains png/pdf of the manuscript figures. 

### **`R`**

R scripts used to generate the objects in `data_processed` and the figure files in `figures`. 

* **`figure-NAME.R`** generates the corresponding figure files `NAME.pdf` or `NAME.png`. The only exception is `figure-roc.R` which generates `roc.pdf` and `prc.pdf`.
* **`helper-functions.R`** small function to calculate peptide ratios.
* **`run-slopes.R`** calculates average antibody response over time across individuals. 
* **`run-ktsp.R`** partitions the data into training/testing sets to find optimal cutoffs and the number of peptide pairs. 
* **`run-validation.R`** evaluates the kTSP classifier on a new data set.
* **`run-cv.R`** uses a cross-validation approach to identifying the number of peptide pairs. 
* **`run-full-data-cutoffs.R`** examines the performance of the classifier if the full data was used to derive the optimal cutoffs for the top performing pairs. 
* **`table-data-summary.R`** generates tables:
    + Table 1: description of the training, testing, and validation sets. 
    + Table S2: description of repeated samples and samples from new visits from individuals in the train+test set.
* **`table-full-data-cutoffs.R`** generates tables:
    + Table S5: comparison of cutoffs identified using the training data only versus the full (train and test) data. 
* **`table-ktsp-results.R`** generates tables:
    + Table 2: performance of the kTSP classifier versus the current LAg + viral load approach. 
    + Table S4: performance of the 3TSP versus 4TSP classifier.
* **`table-pep-summary.R`** generates tables:
    + Table S1: distribution of peptides in spanning the HIV genome.
    

# Mapping of manuscript content to repository files

## Methods

The code for the methods described were run in the order of:

* `R/run-slopes.R`
* `R/run-ktsp.R`
* `R/run-validation.R`
* `R/run-cv.R`
* `R/run-full-data-cutoffs.R`

## Mapping of Figures and Tables

### Main Text

* Table 1 - generated with `R/table-data-summary.R`
* Figure 1 - in `figures/roc.pdf`; generated with `R/figure-roc.R`
* Table 2 - generated with `R/table-ktsp-summary.R`
* Figure 2 - in `figures/logistic.pdf`; generated with `R/figure-logistic.R`

### Supplementary Material

* Figure 1 - in `figures/logfc.png`; generated with `R/figure-logfc.R`
* Figure 2 - in `figures/train-test-k.pdf`; generated with `R/figure-train-test-k.R`
* Figure 3 - in `figures/loocv-k.pdf`; generated with `R/figure-loocv-k.R`
* Figure 4 - in `figures/ktsp_scheme.png`; generated in Microsoft PowerPoint. 
* Figure 5 - in `figures/prc.pdf`; generated with `R/figure-roc.R`
* Table 1 - generated with `R/table-pep-summary.R`
* Table 2 - generated with `R/table-data-summary.R`
* Table 3 - generated with `R/table-pep-summary.R`
* Table 4 - generated with `R/table-ktsp-results.R`
* Table 5 - generated with `R/table-full-data-cutoffs.R`

# Session Info

```
> devtools::session_info()
─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.2 (2020-06-22)
 os       macOS  10.16                
 system   x86_64, darwin17.0          
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       America/New_York            
 date     2020-12-22                  

─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────
 package      * version date       lib source        
 assertthat     0.2.1   2019-03-21 [1] CRAN (R 4.0.0)
 backports      1.2.0   2020-11-02 [1] CRAN (R 4.0.2)
 broom          0.7.2   2020-10-20 [1] CRAN (R 4.0.2)
 callr          3.5.1   2020-10-13 [1] CRAN (R 4.0.2)
 cellranger     1.1.0   2016-07-27 [1] CRAN (R 4.0.0)
 cli            2.2.0   2020-11-20 [1] CRAN (R 4.0.2)
 colorspace     2.0-0   2020-11-11 [1] CRAN (R 4.0.2)
 crayon         1.3.4   2017-09-16 [1] CRAN (R 4.0.0)
 DBI            1.1.0   2019-12-15 [1] CRAN (R 4.0.0)
 dbplyr         2.0.0   2020-11-03 [1] CRAN (R 4.0.2)
 desc           1.2.0   2018-05-01 [1] CRAN (R 4.0.0)
 devtools     * 2.3.2   2020-09-18 [1] CRAN (R 4.0.2)
 digest         0.6.27  2020-10-24 [1] CRAN (R 4.0.2)
 dplyr        * 1.0.2   2020-08-18 [1] CRAN (R 4.0.2)
 ellipsis       0.3.1   2020-05-15 [1] CRAN (R 4.0.0)
 fansi          0.4.1   2020-01-08 [1] CRAN (R 4.0.0)
 farver         2.0.3   2020-01-16 [1] CRAN (R 4.0.0)
 forcats      * 0.5.0   2020-03-01 [1] CRAN (R 4.0.0)
 Formula        1.2-4   2020-10-16 [1] CRAN (R 4.0.2)
 fs             1.5.0   2020-07-31 [1] CRAN (R 4.0.2)
 generics       0.1.0   2020-10-31 [1] CRAN (R 4.0.2)
 ggplot2      * 3.3.2   2020-06-19 [1] CRAN (R 4.0.2)
 glue           1.4.2   2020-08-27 [1] CRAN (R 4.0.2)
 gtable         0.3.0   2019-03-25 [1] CRAN (R 4.0.0)
 haven          2.3.1   2020-06-01 [1] CRAN (R 4.0.0)
 here         * 0.1     2017-05-28 [1] CRAN (R 4.0.2)
 hms            0.5.3   2020-01-08 [1] CRAN (R 4.0.0)
 httr           1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
 jsonlite       1.7.1   2020-09-07 [1] CRAN (R 4.0.2)
 knitr          1.30    2020-09-22 [1] CRAN (R 4.0.2)
 labeling       0.4.2   2020-10-20 [1] CRAN (R 4.0.2)
 lattice        0.20-41 2020-04-02 [1] CRAN (R 4.0.2)
 lifecycle      0.2.0   2020-03-06 [1] CRAN (R 4.0.0)
 lubridate      1.7.9   2020-06-08 [1] CRAN (R 4.0.2)
 magrittr       2.0.1   2020-11-17 [1] CRAN (R 4.0.2)
 Matrix         1.2-18  2019-11-27 [1] CRAN (R 4.0.2)
 memoise        1.1.0   2017-04-21 [1] CRAN (R 4.0.2)
 mgcv           1.8-33  2020-08-27 [1] CRAN (R 4.0.2)
 modelr         0.1.8   2020-05-19 [1] CRAN (R 4.0.0)
 munsell        0.5.0   2018-06-12 [1] CRAN (R 4.0.0)
 nlme           3.1-150 2020-10-24 [1] CRAN (R 4.0.2)
 pillar         1.4.6   2020-07-10 [1] CRAN (R 4.0.2)
 pkgbuild       1.1.0   2020-07-13 [1] CRAN (R 4.0.2)
 pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.0.0)
 pkgload        1.1.0   2020-05-29 [1] CRAN (R 4.0.0)
 plyr           1.8.6   2020-03-03 [1] CRAN (R 4.0.0)
 prettyunits    1.1.1   2020-01-24 [1] CRAN (R 4.0.0)
 processx       3.4.4   2020-09-03 [1] CRAN (R 4.0.2)
 ps             1.4.0   2020-10-07 [1] CRAN (R 4.0.2)
 purrr        * 0.3.4   2020-04-17 [1] CRAN (R 4.0.0)
 R6             2.5.0   2020-10-28 [1] CRAN (R 4.0.2)
 RColorBrewer * 1.1-2   2014-12-07 [1] CRAN (R 4.0.0)
 Rcpp           1.0.5.2 2020-08-13 [1] local         
 readr        * 1.4.0   2020-10-05 [1] CRAN (R 4.0.2)
 readxl         1.3.1   2019-03-13 [1] CRAN (R 4.0.0)
 remotes        2.2.0   2020-07-21 [1] CRAN (R 4.0.2)
 reprex         0.3.0   2019-05-16 [1] CRAN (R 4.0.0)
 rlang          0.4.8   2020-10-08 [1] CRAN (R 4.0.2)
 rprojroot      1.3-2   2018-01-03 [1] CRAN (R 4.0.0)
 rstudioapi     0.13    2020-11-12 [1] CRAN (R 4.0.2)
 rvest          0.3.6   2020-07-25 [1] CRAN (R 4.0.2)
 scales         1.1.1   2020-05-11 [1] CRAN (R 4.0.0)
 sessioninfo    1.1.1   2018-11-05 [1] CRAN (R 4.0.2)
 stringi        1.5.3   2020-09-09 [1] CRAN (R 4.0.2)
 stringr      * 1.4.0   2019-02-10 [1] CRAN (R 4.0.0)
 table1       * 1.2     2020-03-23 [1] CRAN (R 4.0.2)
 testthat       3.0.0   2020-10-31 [1] CRAN (R 4.0.2)
 tibble       * 3.0.4   2020-10-12 [1] CRAN (R 4.0.2)
 tidyr        * 1.1.2   2020-08-27 [1] CRAN (R 4.0.2)
 tidyselect     1.1.0   2020-05-11 [1] CRAN (R 4.0.0)
 tidyverse    * 1.3.0   2019-11-21 [1] CRAN (R 4.0.2)
 usethis      * 1.6.3   2020-09-17 [1] CRAN (R 4.0.2)
 vctrs          0.3.4   2020-08-29 [1] CRAN (R 4.0.2)
 withr          2.3.0   2020-09-22 [1] CRAN (R 4.0.2)
 xfun           0.19    2020-10-30 [1] CRAN (R 4.0.2)
 xml2           1.3.2   2020-04-23 [1] CRAN (R 4.0.0)

[1] /Library/Frameworks/R.framework/Versions/4.0/Resources/library
```