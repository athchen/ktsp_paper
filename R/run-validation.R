# run-validation.R ------------
# Code to look at the ktsp performance on the validation data
# -----------------------

# load discovery data to show duplicates
load("data/hiv_data.rda")
discovery_rc <- sample_anno[, c("key", "ptid", "yrs_post_sero")] %>%
  mutate(key = paste(ptid, yrs_post_sero, sep = "_"))
rm(list = (ls()[ls()!= "discovery_rc"]))

load("data/validation_data.rda")

# Load ktsp optimal cutoffs and peptide pairs
kpep_list <- readRDS("data/kpep_list.rds")

validation_tsp <- rc_ptid_yrs[, c("ptid", "yrs_post_sero", kpep_list$pep_1, kpep_list$pep_2)] %>%
  mutate(`pep_25409-pep_77241` = pep_25409 - pep_77241, 
         `pep_93864-pep_92687` = pep_93864 - pep_92687, 
         `pep_26777-pep_20904` = pep_26777 - pep_20904, 
         `pep_77323-pep_32119` = pep_77323 - pep_32119) %>%
  dplyr::select(ptid, yrs_post_sero, kpep_list$pep_pair[1:4]) %>%
  gather("pep_pair", "rc", -ptid, -yrs_post_sero) %>%
  left_join(kpep_list[, c("pep_pair", "opti_cutoff")], by = "pep_pair") %>%
  mutate(ktsp_votes = ifelse(rc <= opti_cutoff, 1, 0)) %>%
  dplyr::select(ptid, yrs_post_sero, pep_pair, ktsp_votes) %>%
  group_by(ptid, yrs_post_sero) %>%
  summarize(ktsp_rec = ifelse(sum(ktsp_votes) >= 3, 1, 0), .groups = "drop") %>%
  ungroup() %>%
  mutate(acc_rec = ifelse(yrs_post_sero <= 1, 1, 0))

# saveRDS(validation_tsp, file = "data_processed/validation_tsp.rds")