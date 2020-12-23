# table-pep-summary.R ----------
# Code to generate:
#   * Table 1 in the supplement
#   * Table 3 in the supplement
# -------------------------------

load(here("data_raw/hiv_data.rda"))

# Table S1 -----------------------
pep_anno %>%
  filter(pep_id %in% pep_anno$pep_id[hiv_ind]) %>%
  group_by(species) %>%
  summarize(num_peptides = n(), .groups = "drop") %>%
  mutate(species = gsub("Primate lentivirus group Human immunodeficiency virus type ", "HIV-", species))

# Table S3 -----------------------
kpep_list <- readRDS(here("data_processed/kpep_list.rds"))

kpep_list %>% 
  mutate(pair = 1:10) %>%
  pivot_longer(cols = c("pep_1", "pep_2"), 
               names_to = "pep_lab", 
               values_to = "peptide") %>%
  mutate(pep_lab = ifelse(pep_lab == "pep_1", "A", "B")) %>%
  unite("peptide_lab", c("pair", "pep_lab"), sep = "") %>%
  left_join(pep_anno, by = c("peptide" = "pep_id")) %>%
  mutate(subtype = str_match(species, "subtype (.*)")[, 2]) %>%
  select(peptide_lab, subtype, pep_pos, pep_aa, uniprot_acc, opti_cutoff)