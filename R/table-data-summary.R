# table-data-summary.R ----------
# Code to generate:
#   * Table 1 in the main text
#   * Table 2 in the supplement
# -------------------------------

library(table1)

load(here("data_raw/hiv_data.rda"))

rndr.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=4), 
       c("", "Median [Min, Max]"=sprintf("%s [%s, %s]", MEDIAN, MIN, MAX)))
}

# Table 1 -----------------------
## Training and testing data
RNGversion("3.5.3")
set.seed(369)
ptid_train  <- sample(pt_anno$ptid, 38, replace = FALSE)
ptid_test   <- pt_anno$ptid[!(pt_anno$ptid %in% ptid_train)]

discovery <- sample_anno %>%
  filter((yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5) &
           yrs_post_sero >= 1/6) %>%
  mutate(set = factor(ifelse(ptid %in% ptid_train, "train", "test"), 
                      levels = c("train", "test")), 
         recent = factor(ifelse(yrs_post_sero <= 1, "recent", "non-recent"), 
                         levels = c("recent", "non-recent"))) %>%
  dplyr::select(set, ptid, recent, yrs_post_sero, age_at_sero, cd4, viral_load) 

label(discovery$recent) <- "Number of samples"
label(discovery$yrs_post_sero) <- "Duration of infection"
label(discovery$age_at_sero) <- "Age of participant"
label(discovery$cd4) <- "CD4 cell count (cells/mm3)"
label(discovery$viral_load) <- "log10 viral load (copies/mL)"

# per person summary
discovery %>%
  group_by(ptid, set) %>%
  summarize(num_samples = n(), .groups = "drop") %>%
  group_by(set) %>%
  summarize(num_people = n(), 
            samples_median = median(num_samples), 
            sample_min = min(num_samples), 
            sample_max = max(num_samples)) 

# duration of infection, age, cd4, viral load summary
table1::table1(~ yrs_post_sero + age_at_sero + cd4 + viral_load | recent*set, 
               discovery, render.continuous = rndr.cont)

## Validation data
rm(list = ls()[!ls() %in% c("discovery", "rndr.cont", "ptid_train", "ptid_test")])
load(here("data_raw/validation_data.rda"))

validation <- sample_anno %>%
  left_join(pt_anno %>% select(ptid, subtype), by = c("ptid")) %>%
  filter((yrs_post_sero <= 0.5 | yrs_post_sero >= 1.5) &
           yrs_post_sero >= 1/6 &
           subtype == "C") %>%
  mutate(set = case_when(paste(ptid, round(yrs_post_sero, 2), sep = "_") %in%
                           paste(discovery$ptid, discovery$yrs_post_sero, sep = "_") ~ 
                           "Repeated individuals, replicated samples", 
                         ptid %in% discovery$ptid ~ 
                           "Repeated individuals, different study visits",
                         TRUE ~ "Validation set"), 
         recent = factor(ifelse(yrs_post_sero <= 1, "recent", "non-recent"), 
                         levels = c("recent", "non-recent"))) %>%
  mutate(set = factor(set, levels = c("Validation set", 
                                      "Repeated individuals, replicated samples",
                                      "Repeated individuals, different study visits")))

label(validation$recent) <- "Number of samples"
label(validation$yrs_post_sero) <- "Duration of infection"
label(validation$age_at_sero) <- "Age of participant"
label(validation$cd4) <- "CD4 cell count (cells/mm3)"
label(validation$viral_load) <- "log10 viral load (copies/mL)"

# per person summary
validation %>%
  group_by(ptid, set) %>%
  summarize(num_samples = n(), .groups = "drop") %>%
  group_by(set) %>%
  summarize(num_people = n(), 
            samples_median = median(num_samples), 
            sample_min = min(num_samples), 
            sample_max = max(num_samples)) 

# duration of infection, age, cd4, viral load summary
table1::table1(~ yrs_post_sero + age_at_visit + cd4 + viral_load | recent, 
               validation %>% filter(set == "Validation set"), render.continuous = rndr.cont)

# Table S2 ----------------------
table1::table1(~ yrs_post_sero + age_at_visit + cd4 + viral_load | recent*set, 
               validation %>% filter(set != "Validation set"), render.continuous = rndr.cont)
