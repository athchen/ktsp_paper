# figure-logfc.R ----------
# Code to generate:
#   * `figures/logfc.png`
# -----------------------

load(here("data_raw/hiv_data.rda"))
kpep_list <- readRDS(here("data_processed/kpep_list.rds"))

rc_ptid_yrs[, c("ptid", "yrs_post_sero", kpep_list$pep_1[1:2], kpep_list$pep_2[1:2])] %>%
  filter(yrs_post_sero >= 1/6) %>%
  pivot_longer(cols = contains("pep"), 
               names_to = "peptide", 
               values_to = "log10_rc") %>%
  left_join(data.frame(peptide = c(kpep_list$pep_1[1:2], kpep_list$pep_2[1:2]), 
                       pep_alias = c("Peptide 1A", "Peptide 2A", "Peptide 1B", "Peptide 2B")), 
            by = c("peptide")) %>%
  ggplot(aes(x = log10(yrs_post_sero), y = log10_rc)) +
  geom_line(aes(group = ptid), color = "grey", size = 0.2) + 
  facet_wrap(pep_alias ~., nrow = 2) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "duration of infection", 
       y = "log10 normalized signal") +
  scale_x_continuous(breaks = log10(c(1/6, 1/2, 1, 2, 5)), 
                     labels = c("2mo", "6mo", "1yr", "2yr", "5yrs")) +
  theme_bw()

ggsave("figures/logfc.png", units = "in", width = 6, height = 6)
  