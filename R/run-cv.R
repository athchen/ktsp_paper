# run-cv.R ----------
# Code to run leave-one-out cross-validation (LOOCV) 
# to select k
# -----------------------

load(here("data_raw/hiv_data.rda"))
kpep_list <- readRDS(here("data_processed/kpep_list.rds"))
source(here("R/helper-functions.R"))

# Make df of ratios for all 10 peptide pairs as well as the predictions
pep_ratios <- apply(kpep_list, 1, function(row){
  get_ratio(row[["pep_1"]], row[["pep_2"]], rc_ptid_yrs)
}) %>% plyr::ldply() %>%
  left_join(kpep_list %>%
              mutate(pair_num = row_number()) %>%
              dplyr::select(pair_num, pep_pair, opti_cutoff), 
            by = c("pep_pair")) %>%
  rowwise() %>%
  mutate(pred_recent = ifelse(logratio <= opti_cutoff, 1, 0), 
         true_recent = ifelse(yrs_post_sero <= 1, 1, 0)) %>%
  filter((yrs_post_sero >= 1/6 & yrs_post_sero <= 0.5) |
           (yrs_post_sero >= 1.5)) 

# make the kTSP prediction matrix using k = 1, 2, ..., 10 pairs. 
ktsp_votes <- pep_ratios %>%
  group_by(ptid, yrs_post_sero, true_recent) %>%
  group_split() %>%
  map_df(function(df){
    num_votes <- sapply(1:10, function(k){
      sum(df$pred_recent[df$pep_pair %in% kpep_list$pep_pair[1:k]])
    })
    
    data.frame(ptid = df$ptid, 
               yrs_post_sero = df$yrs_post_sero, 
               true_recent = df$true_recent, 
               num_pairs = 1:10,
               num_votes = num_votes) %>%
      mutate(kpred = ifelse(num_votes > num_pairs/2, 1, 0))
  })


# ---- Cross-validation -----
# get cv partitions
cv_list <- data.frame(group = 1:57, ptids = pt_anno$ptid)

# For each partition of the ptids, we do the following:
# The training accuracies for all pairs along with the final outcome 
# are returned as a list with 2 elements. 
cv_results <- lapply(1:57, function(test_set){
  # get ptids in training set
  ptid_train <- cv_list$ptids[cv_list$group != test_set]
  ptid_test <- cv_list$ptids[cv_list$group == test_set]
  
  # get classification accuracy for k = 1, 2, ..., 10
  accuracy_df <- ktsp_votes %>% 
    filter(ptid %in% ptid_train) %>%
    group_by(num_pairs) %>%
    summarize(accuracy = sum(true_recent == kpred)/n(), 
              .groups = "drop")
  
  # get optimal number of pairs, pick the smallest number if two are present. 
  opti_pairs <- accuracy_df %>%
    filter(accuracy == max(accuracy)) %>%
    pull(num_pairs) %>%
    as.numeric() 
  
  # calculate test accuracy
  test_accuracy <- ktsp_votes %>%
    filter(ptid %in% ptid_test) %>%
    group_by(num_pairs) %>%
    summarize(accuracy = sum(true_recent == kpred)/n(), 
              num_samples = n(), 
              .groups = "drop") 
  
  return(list(train_results = accuracy_df %>%
                mutate(test_sample = ptid_test), 
              opti_pairs = tibble(opti_pairs = opti_pairs, 
                                  test_sample = ptid_test), 
              test_results  = test_accuracy %>%
                mutate(test_sample = ptid_test)))
})

# ----- tidy results -----
cv_train <- cv_results %>%
  map_dfr(function(y){y$train_results}) %>%
  mutate(num_pairs = as.integer(num_pairs))

cv_opti <- cv_results %>%
  map_dfr(function(y){y$opti_pairs})

cv_test <- cv_results %>%
  map_dfr(function(y){y$test_results}) %>%
  mutate(num_pairs = as.integer(num_pairs))
