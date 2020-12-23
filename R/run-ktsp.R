# run-ktsp.R ------------
# Code to form peptide pairs, identify optimal cutoffs, and select k 
# using a train/test partition. 
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

# Calculate all peptide ratios
pep_ratios <- apply(pairs, 1, function(row) get_ratio(row[1], row[2], rc_ptid_yrs)) %>%
  plyr::ldply(.id = NULL) %>%
  spread(pep_pair, logratio) %>%
  arrange(ptid)

# Split into training and testing set -----------------------
RNGversion("3.5.3")
set.seed(369)
ptid_train  <- sample(pt_anno$ptid, 38, replace = FALSE)
ptid_test   <- pt_anno$ptid[!(pt_anno$ptid %in% ptid_train)]

# Split data into train and test, and dichotomize into recent/non-recent
#   Recent = 1
#   Non-recent = 0
pep_ratio_train <- pep_ratios %>%
  filter(ptid %in% ptid_train &
           yrs_post_sero >= 1/6 & 
           (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5)) %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0))

pep_ratio_test <- pep_ratios %>%
  filter(ptid %in% ptid_test &
           yrs_post_sero >= 1/6 &
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


training_cutoffs <- lapply(rownames(pairs), function(x) get_cutoffs(pep_ratio_train, x)) 

# data for ROC/PRC curves
train_curve <- training_cutoffs %>%
  lapply(function(x) x[[1]]) %>%
  plyr::ldply()

# df of optimal cutofs for each peptide pairs
train_opti_cutoffs <- training_cutoffs %>%
  lapply(function(x) x[[2]]) %>%
  plyr::ldply()

# df of pep_pair votes based on optimal cutoffs
train_votes <- pep_ratio_train %>%
  pivot_longer(cols = contains("pep"),
               names_to = "pep_pair", 
               values_to = "rc") %>%
  left_join(train_opti_cutoffs %>%
              select(pep_pair, opti_cutoff),
            by = "pep_pair") %>%
  mutate(pred_rec = ifelse(rc<= opti_cutoff, 1, 0)) 

# Order peptide pairs -----------------------
# Using the optimal cutoffs for each peptide pair found above, 
# create a set of peptide pairs as follows: 
#   1. Take the peptide pair with the higest accuracy. 
#       If there are any ties, we pick the peptide pair 
#       with the largest area under the precision-recall curve. 
#   2. Remove any peptide pairs that share a peptide with the 
#       peptide pairs in the list. Repeat step 1 until there 
#       are no remaining peptide pairs. 

kpep_list <- data.frame()
pairs_df <- pairs %>% rownames_to_column("pep_pair") %>%
  rename(pep_1 = Var1, 
         pep_2 = Var2)

train_accuracy <- train_votes %>%
  group_by(pep_pair) %>%
  summarize(perc_correct = sum((recent == 1 & pred_rec == 1) |
                                 (recent == 0 & pred_rec == 0))/length(pep_pair), 
            window_perc = sum(((recent == 1 & pred_rec == 1) |
                                 (recent == 0 & pred_rec == 0)) &
                                (yrs_post_sero >= 1.5 | yrs_post_sero <= 0.5))/
              sum(yrs_post_sero >= 1.5 | yrs_post_sero <= 0.5), 
            .groups = "drop") %>%
  left_join(train_opti_cutoffs %>% 
              dplyr::select(pep_pair, opti_cutoff), by = "pep_pair") %>%
  left_join(train_curve %>% 
              dplyr::select(pep_pair, auc_prc) %>%
              distinct(), by = "pep_pair")

while(dim(pairs_df)[1] > 0){
  
  # Subset data to remaining pairs
  data_sub <- train_accuracy %>% filter(pep_pair %in% pairs_df$pep_pair)
  
  # Select top peptide pair
  top_pairs <- data_sub %>% filter(window_perc == max(window_perc))
  
  # Check if there are more than one peptide pair
  if(dim(top_pairs)[1] != 1) {
    best_pair <- top_pairs %>% filter(auc_prc == max(auc_prc))
  } else{
    best_pair <- top_pairs
  }
  
  # Check that there is only one pair left
  if(dim(best_pair)[1] != 1){
    message("More than one pair with the same auc_prc")
    break
  }
  
  # Add best pair to the list
  best_pair_names <- best_pair %>% 
    separate(pep_pair, c("pep_1", "pep_2"), sep = "-", remove = F) %>%
    dplyr::select(pep_pair, pep_1, pep_2) 
  
  kpep_list <- rbind(kpep_list, best_pair_names)
  
  # Remove peptides that are common with the peptide pair
  pairs_df <- pairs_df %>% filter(!(pep_1 %in% kpep_list$pep_1) &
                                        !(pep_2 %in% kpep_list$pep_2))
}

# kpep_list %>%
#   left_join(train_opti_cutoffs, by = c("pep_pair")) %>%
#   select(pep_pair, pep_1, pep_2, opti_cutoff) %>%
#   saveRDS(., file = "data/kpep_list.rds")

# Select k based on training and testing data -----------------------
train_ktsp <- train_votes %>%
  filter(pep_pair %in% kpep_list$pep_pair) %>%
  group_by(ptid, yrs_post_sero, recent) %>%
  group_split() %>%
  map_df(function(df){
    num_votes <- sapply(1:10, function(k){
      sum(df$pred_rec[df$pep_pair %in% kpep_list$pep_pair[1:k]])
    })
    
    data.frame(ptid = df$ptid, 
               yrs_post_sero = df$yrs_post_sero, 
               recent = df$recent, 
               num_pairs = 1:10,
               num_votes = num_votes) %>%
      mutate(kpred = ifelse(num_votes > num_pairs/2, 1, 0))
  }) 

# Testing
test_ktsp <- pep_ratio_test %>%
  pivot_longer(cols = contains("pep"), 
               names_to = "pep_pair", 
               values_to = "rc") %>%
  left_join(train_opti_cutoffs, by = c("pep_pair")) %>%
  mutate(pred_rec = ifelse(rc <= opti_cutoff, 1, 0)) %>%
  filter(pep_pair %in% kpep_list$pep_pair) %>%
  group_by(ptid, yrs_post_sero, recent) %>%
  group_split() %>%
  map_df(function(df){
    num_votes <- sapply(1:10, function(k){
      sum(df$pred_rec[df$pep_pair %in% kpep_list$pep_pair[1:k]])
    })
    
    data.frame(ptid = df$ptid, 
               yrs_post_sero = df$yrs_post_sero, 
               recent = df$recent, 
               num_pairs = 1:10,
               num_votes = num_votes) %>%
      mutate(kpred = ifelse(num_votes > num_pairs/2, 1, 0))
  }) 

# Get LAg performance -----------------------
# Create data frame of only window samples
lag_train <- sample_anno %>%
  filter(yrs_post_sero >= 1/6 &
           (ptid %in% ptid_train) &
           (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5)) %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0)) %>%
  dplyr::select(ptid, yrs_post_sero, recent, lag)

# Get lag curve
lag_curve <- get_cutoffs(lag_train, "lag")