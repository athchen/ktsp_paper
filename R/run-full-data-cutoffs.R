# run-full-data-cutoffs.R ------------
# Code to identify peptide cutoffs using the full data (train + test)
# -----------------------

load(here("data_raw/hiv_data.rda"))
load(here("data_processed/slopes.rda"))
source(here("R/helper-functions.R"))

# Identify top 10 increasing peptides
top10_upgoing <- slopes %>%
  group_by(pep_id) %>%
  summarize(avg_slope = mean(slope), 
            .groups = "drop") %>%
  top_n(10, avg_slope) %>%
  pull(pep_id) %>%
  as.character()

# Identify top 10 decreasing peptides
top10_downgoing <- slopes %>%
  group_by(pep_id) %>%
  summarize(avg_slope = mean(slope), 
            .groups = "drop") %>%
  top_n(-10, avg_slope) %>%
  pull(pep_id) %>%
  as.character()

# Create matrix of peptide pairs
pairs <- expand.grid(top10_upgoing, top10_downgoing)
rownames(pairs) <- paste0(pairs[,1], "-", pairs[, 2])

# Calculate all peptide ratios for training + testing data
pep_ratios <- apply(pairs, 1, function(row) get_ratio(row[1], row[2], rc_ptid_yrs)) %>%
  plyr::ldply(.id = NULL) %>%
  spread(pep_pair, logratio) %>%
  arrange(ptid) %>%
  filter(yrs_post_sero >= 1/6 & 
           (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5)) %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0))

# Get optimal cutoffs -----------------------
#' Function to calculate the ppv, sens, spec, AUC for ROC and PRC, 
#' and identify optimal cutoffs for each peptide pair
#'
#' @param data data frame with `ptid`, `yrs_post_sero`, and the peptides in columns
#' @param peptide peptide pair of interest
#' @param extra_info boolean indicating whether there is a column with classification.
#' The column must be named `extra_info`
#'
#' @return list of three objecs:
#' * `pep_all`: data frame with ppvs, sens, spec, auc_roc, and suc_prc for each cutoff
#' * `pep_sum`: data frame with the peptide pair, optimal cutoff, # of cutoffs with the
#'  same optimality conditions, and the min/max/range of the optimal cutoffs
#' * `cutoffs`: list of cutoffs with the same ranking for each peptide pair
get_cutoffs <- function(data, peptide, extra_info = FALSE){
  
  # Find range of peptide RC, then create cutoffs by partitioning the range into 1000 points
  min_cutoff <- min(data[, peptide])
  max_cutoff <- max(data[, peptide])
  cutoffs <- seq(min_cutoff, max_cutoff, length.out = 1000)
  
  # Calculate ppv, sens, spec for each cutoff
  ppv <- sapply(cutoffs, function(x) {
    if(extra_info){
      predict <- (data[, peptide] <= x & data$extra_info)
    } else {
      predict <- (data[, peptide] <= x)
    }
    
    sum((data$recent == 1) & predict)/sum(predict)
  })
  
  spec <- sapply(cutoffs, function(x){
    if(extra_info){
      predict <- (data[, peptide] <= x & data$extra_info)
    } else {
      predict <- (data[, peptide] <= x)
    }
    
    sum((data$recent == 0) & !predict)/sum(data$recent == 0)
  })
  
  sens <- sapply(cutoffs, function(x){
    if(extra_info){
      predict <- (data[, peptide] <= x & data$extra_info)
    } else {
      predict <- (data[, peptide] <= x)
    }
    sum((data$recent == 1) & predict)/sum(data$recent == 1)
  })
  
  # Get AUC for ROC
  tmp <- data.frame(sens = sens, spec_c = 1-spec) %>%
    arrange(spec_c)
  npoints <- dim(tmp)[1]
  area <- sum(0.5 * (tmp$sens[-1] + tmp$sens[-npoints]) * (tmp$spec_c[-1] -
                                                             tmp$spec_c[-npoints]))
  
  # Get AUC for PRC
  tmp <- data.frame(sens = sens, ppv = ppv) %>%
    arrange(sens)
  npoint <- dim(tmp)[1]
  area_prc <- sum(0.5 * (tmp$ppv[-1] + tmp$ppv[-npoints]) * (tmp$sens[-1] - tmp$sens[-npoints]))
  
  # Get optimal cutoff using the minimum distance between PRC curve and (1,1),
  # Take the middle cutoff. If the number of cutoffs is even,
  # it takes the larger of the two middle cutoffs.
  d <- sqrt((1-sens)^2 + (1-ppv)^2)
  # d <- sqrt((0-(1-spec))^2 + (1 - sens)^2)
  opti_cutoff <- cutoffs[d == min(d)][ceiling(length(cutoffs[d == min(d)])/2)]
  
  return(list(pep_all = data.frame(pep_pair = peptide,
                                   cutoff = cutoffs,
                                   ppv = ppv,
                                   sens = sens,
                                   spec = spec,
                                   auc_roc = area,
                                   auc_prc = area_prc), 
              pep_sum = data.frame(pep_pair = peptide,
                                   opti_cutoff = opti_cutoff, 
                                   n_cutoffs = sum(d == min(d)), 
                                   min_cutoffs = min(cutoffs[d == min(d)]), 
                                   max_cutoffs = max(cutoffs[d == min(d)]), 
                                   range = max(cutoffs[d == min(d)]) - min(cutoffs[d == min(d)])), 
              cutoffs = cutoffs[d == min(d)]))
}

cutoffs <- lapply(rownames(pairs), function(x) get_cutoffs(pep_ratios, x)[[2]]) %>%
  plyr::ldply()
