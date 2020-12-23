# table-ktsp-summary.R ----------
# Code to generate:
#   * Table 2 in the main text
#   * Table 4 in the supplement
# -------------------------------

# load validation results
source(here("R/run-validation.R"))
validation_anno <- sample_anno %>%
  mutate(lag_prediction = ifelse(lag >= 1.5 | (lag < 1.5 & viral_load <= 3), 0, 1)) %>%
  left_join(pt_anno %>% select(ptid, subtype), by = c("ptid"))
  
# load train/test results
source(here("R/run-ktsp.R"))
discovery_anno <- sample_anno %>%
  mutate(lag_prediction = ifelse(lag >= 1.5 | (lag < 1.5 & viral_load <= 3), 0, 1))

rm(list = ls()[!ls() %in% c("kpep_list", "discovery_anno", "validation_anno", 
                            "train_ktsp", "test_ktsp", "validation_tsp",
                            "ptid_train", "ptid_test")])

# Table 2 -----------------------
lag_merged <- discovery_anno %>%
  filter(yrs_post_sero >= 1/6 &
           (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5)) %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0), 
         data = ifelse(ptid %in% ptid_train, "train", "test")) %>%
  select(ptid, yrs_post_sero, lag_prediction, data) %>%
  bind_rows(validation_anno %>%
              filter(yrs_post_sero >= 1/6 &
                       (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5) &
                       subtype == "C") %>%
              mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0), 
                     data = case_when(paste(ptid, round(yrs_post_sero, 2), sep = "_") %in%
                                        paste(discovery_anno$ptid, discovery_anno$yrs_post_sero, sep = "_") ~ 
                                        "Repeated individuals, replicated samples", 
                                      ptid %in% discovery_anno$ptid ~ 
                                        "Repeated individuals, different study visits",
                                      TRUE ~ "Validation set"))) %>%
              select(ptid, yrs_post_sero, lag_prediction, data)

table2_dat <- bind_rows(train_ktsp %>% 
                          filter(num_pairs == 4) %>%
                          mutate(data = "train"), 
                        test_ktsp %>% 
                          filter(num_pairs == 4) %>%
                          mutate(data = "test")) %>%
  select(ptid, yrs_post_sero, recent, kpred, data) %>%
  bind_rows(validation_tsp %>%
              left_join(validation_anno %>% select(ptid, subtype) %>% distinct(),
                        by = c("ptid")) %>%
              filter(subtype == "C" &
                       yrs_post_sero >= 1/6 &
                       (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5)) %>%
              rename(recent = acc_rec, 
                     kpred = ktsp4_rec) %>%
              mutate(data = case_when(paste(ptid, round(yrs_post_sero, 2), sep = "_") %in%
                                        paste(discovery_anno$ptid, discovery_anno$yrs_post_sero, sep = "_") ~ 
                                        "Repeated individuals, replicated samples", 
                                      ptid %in% discovery_anno$ptid ~ 
                                        "Repeated individuals, different study visits",
                                      TRUE ~ "Validation set")) %>%
              select(-subtype, -ktsp3_rec)) %>%
  left_join(lag_merged, by = c("ptid", "yrs_post_sero", "data"))

# Totals
table2_dat %>%
  pivot_longer(cols = c("kpred", "lag_prediction"), 
               names_to = "method", 
               values_to = "prediction") %>%
  group_by(data, method) %>%
  summarize(correct = sum(recent == prediction, na.rm = TRUE), 
            total = sum(!is.na(prediction)), 
            percent = correct/total, .groups = "drop") %>%
  mutate(data = factor(data, levels = c("train", "test", "Validation set", 
                                        "Repeated individuals, replicated samples", 
                                        "Repeated individuals, different study visits"))) %>%
  arrange(data) %>%
  pivot_wider(names_from = method, 
              values_from = c("correct", "total", "percent"), 
              names_glue = "{method}_{.value}") %>%
  select(data, kpred_correct, kpred_total, kpred_percent, 
         lag_prediction_correct, lag_prediction_total, lag_prediction_percent)

# Split by recent/non-recent
table2_dat %>%
  pivot_longer(cols = c("kpred", "lag_prediction"), 
               names_to = "method", 
               values_to = "prediction") %>%
  group_by(data, recent, method) %>%
  summarize(correct = sum(recent == prediction, na.rm = TRUE), 
            total = sum(!is.na(prediction)), 
            percent = correct/total, .groups = "drop") %>%
  mutate(recent = factor(ifelse(recent == 1, "recent", "non-recent"), 
                         levels = c("recent", "non-recent")), 
         data = factor(data, levels = c("train", "test", "Validation set", 
                                        "Repeated individuals, replicated samples", 
                                        "Repeated individuals, different study visits"))) %>%
  arrange(data, recent) %>%
  pivot_wider(names_from = method, 
              values_from = c("correct", "total", "percent"), 
              names_glue = "{method}_{.value}") %>%
  select(data, kpred_correct, kpred_total, kpred_percent, 
         lag_prediction_correct, lag_prediction_total, lag_prediction_percent)

# Table S4 ----------------------
table4_dat <- bind_rows(train_ktsp %>% 
                          filter(num_pairs %in% 3:4) %>%
                          mutate(data = "train"), 
                        test_ktsp %>% 
                          filter(num_pairs %in% 3:4) %>%
                          mutate(data = "test")) %>%
  select(ptid, yrs_post_sero, recent, num_pairs, kpred, data) %>%
  bind_rows(validation_tsp %>%
              left_join(validation_anno %>% select(ptid, subtype) %>% distinct(),
                        by = c("ptid")) %>%
              filter(subtype == "C" &
                       yrs_post_sero >= 1/6 &
                       (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5)) %>%
              pivot_longer(cols = c("ktsp3_rec", "ktsp4_rec"), 
                           names_to = "num_pairs", 
                           names_pattern = "ktsp([34])_rec",
                           values_to = "kpred") %>%
              rename(recent = acc_rec) %>%
              mutate(num_pairs = as.numeric(num_pairs), 
                     data = case_when(paste(ptid, round(yrs_post_sero, 2), sep = "_") %in%
                                        paste(discovery_anno$ptid, discovery_anno$yrs_post_sero, sep = "_") ~ 
                                        "Repeated individuals, replicated samples", 
                                      ptid %in% discovery_anno$ptid ~ 
                                        "Repeated individuals, different study visits",
                                      TRUE ~ "Validation set")) %>%
              select(-subtype)) 

# Totals
table4_dat %>%
  group_by(data, num_pairs) %>%
  summarize(correct = sum(recent == kpred, na.rm = TRUE), 
            total = sum(!is.na(kpred)), 
            percent = correct/total, .groups = "drop") %>%
  mutate(data = factor(data, levels = c("train", "test", "Validation set", 
                                        "Repeated individuals, replicated samples", 
                                        "Repeated individuals, different study visits"))) %>%
  arrange(data) %>%
  pivot_wider(names_from = num_pairs, 
              values_from = c("correct", "total", "percent"), 
              names_glue = "pairs{num_pairs}_{.value}") %>%
  select(data, pairs3_correct, pairs3_total, pairs3_percent, 
         pairs4_correct, pairs4_total, pairs4_percent)

# Split by recent/non-recent
table4_dat %>%
  group_by(data, recent, num_pairs) %>%
  summarize(correct = sum(recent == kpred, na.rm = TRUE), 
            total = sum(!is.na(kpred)), 
            percent = correct/total, .groups = "drop") %>%
  mutate(recent = factor(ifelse(recent == 1, "recent", "non-recent"), 
                         levels = c("recent", "non-recent")), 
         data = factor(data, levels = c("train", "test", "Validation set", 
                                        "Repeated individuals, replicated samples", 
                                        "Repeated individuals, different study visits"))) %>%
  arrange(data, recent) %>%
  pivot_wider(names_from = num_pairs, 
              values_from = c("correct", "total", "percent"), 
              names_glue = "pairs{num_pairs}_{.value}") %>%
  select(data, pairs3_correct, pairs3_total, pairs3_percent, 
         pairs4_correct, pairs4_total, pairs4_percent)