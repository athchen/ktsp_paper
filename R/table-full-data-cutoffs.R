# table-full-data-cutoffs.R -----
# Code to generate:
#   * Table 5 in the supplement
# -------------------------------

source(here("R/helper-functions.R"))

kpep_list <- readRDS(here("data_processed/kpep_list.rds"))
source(here("R/run-full-data-cutoffs.R"))
kpep_list <- kpep_list %>%
  left_join(cutoffs %>% 
              rename(full_cutoff = opti_cutoff) %>%
              select(pep_pair, full_cutoff), by = c("pep_pair"))

# Get training/testing ptid
RNGversion("3.5.3")
set.seed(369)
ptid_train  <- sample(pt_anno$ptid, 38, replace = FALSE)
ptid_test   <- pt_anno$ptid[!(pt_anno$ptid %in% ptid_train)]

# Tidy discovery results
pep_ratios <- apply(kpep_list[1:4, ], 1, function(row) get_ratio(row[["pep_1"]], row[["pep_2"]], rc_ptid_yrs)) %>%
  plyr::ldply() 

# discovery outcomes
discovery_outcomes <- pep_ratios %>%
  left_join(kpep_list %>%
              select(pep_pair, opti_cutoff, full_cutoff), by = c("pep_pair")) %>%
  rename(tt_cutoff = opti_cutoff) %>%
  pivot_longer(cols = c("tt_cutoff", "full_cutoff"), 
               names_to = "cutoff_type", 
               values_to = "cutoff", 
               names_pattern = "(.*)_cutoff") %>%
  mutate(pair_vote = ifelse(logratio <= cutoff, 1, 0)) %>%
  group_by(ptid, yrs_post_sero, cutoff_type) %>%
  summarize(ktsp_pred = ifelse(sum(pair_vote) >= 3, 1, 0), .groups = "drop") %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0), 
         data = ifelse(ptid %in% ptid_train, "train", "test")) %>%
  filter(yrs_post_sero >= 1/6 &
           (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5))

discovery_anno <- sample_anno

# validation_outcomes
rm(list = ls()[!ls() %in% c("kpep_list", "get_ratio", "discovery_outcomes", 
                           "discovery_anno")])
load(here("data_raw/validation_data.rda"))
validation_anno <- sample_anno

pep_ratios <- apply(kpep_list[1:4, ], 1, function(row) get_ratio(row[["pep_1"]], row[["pep_2"]], rc_ptid_yrs)) %>%
  plyr::ldply()

validation_outcomes <- pep_ratios %>%
  left_join(kpep_list %>%
              select(pep_pair, opti_cutoff, full_cutoff), by = c("pep_pair")) %>%
  rename(tt_cutoff = opti_cutoff) %>%
  pivot_longer(cols = c("tt_cutoff", "full_cutoff"), 
               names_to = "cutoff_type", 
               values_to = "cutoff", 
               names_pattern = "(.*)_cutoff") %>%
  mutate(pair_vote = ifelse(logratio <= cutoff, 1, 0)) %>%
  group_by(ptid, yrs_post_sero, cutoff_type) %>%
  summarize(ktsp_pred = ifelse(sum(pair_vote) >= 3, 1, 0), .groups = "drop") %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0), 
         data = case_when(paste(ptid, round(yrs_post_sero, 2), sep = "_") %in%
                            paste(discovery_anno$ptid, discovery_anno$yrs_post_sero, sep = "_") ~ 
                            "Repeated individuals, replicated samples", 
                          ptid %in% discovery_anno$ptid ~ 
                            "Repeated individuals, different study visits",
                          TRUE ~ "Validation set")) %>%
  filter(yrs_post_sero >= 1/6 &
           (yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5)) %>%
  left_join(pt_anno %>% select(ptid, subtype), by = "ptid") %>%
  filter(subtype == "C") %>%
  select(-subtype)

# Table S5 ----------------------
# Totals
bind_rows(discovery_outcomes, validation_outcomes) %>%
  group_by(data, cutoff_type) %>%
  summarize(correct = sum(recent == ktsp_pred, na.rm = TRUE), 
            total = sum(!is.na(ktsp_pred)), 
            percent = correct/total, .groups = "drop") %>%
  mutate(data = factor(data, levels = c("train", "test", "Validation set", 
                                        "Repeated individuals, replicated samples", 
                                        "Repeated individuals, different study visits"))) %>%
  arrange(data) %>%
  pivot_wider(names_from = cutoff_type, 
              values_from = c("correct", "total", "percent"), 
              names_glue = "cutoff_{cutoff_type}_{.value}") %>%
  select(data, cutoff_tt_correct, cutoff_tt_total, cutoff_tt_percent,
         cutoff_full_correct, cutoff_full_total, cutoff_full_percent)

# Split by recent/non-recent
bind_rows(discovery_outcomes, validation_outcomes) %>%
  group_by(data, recent, cutoff_type) %>%
  summarize(correct = sum(recent == ktsp_pred, na.rm = TRUE), 
            total = sum(!is.na(ktsp_pred)), 
            percent = correct/total, .groups = "drop") %>%
  mutate(recent = factor(ifelse(recent == 1, "recent", "non-recent"), 
                         levels = c("recent", "non-recent")), 
         data = factor(data, levels = c("train", "test", "Validation set", 
                                        "Repeated individuals, replicated samples", 
                                        "Repeated individuals, different study visits"))) %>%
  arrange(data, recent) %>%
  pivot_wider(names_from = cutoff_type, 
              values_from = c("correct", "total", "percent"), 
              names_glue = "cutoff_{cutoff_type}_{.value}") %>%
  select(data, cutoff_tt_correct, cutoff_tt_total, cutoff_tt_percent,
         cutoff_full_correct, cutoff_full_total, cutoff_full_percent)