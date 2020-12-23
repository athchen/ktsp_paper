# run-slopes.R ----------
# Code to fit a simple linear model of log-normalized relative 
# fold changes as a function of log time for each peptide and each study participant 
# using all available samples collected more than 2 months after infection.
# 
# Generates `data_processed/slopes.rda`
# -----------------------

load(here("data_raw/hiv_data.rda"))

#' Helper function to calculate to tidy and fit the linear model
#'
#' @param pt_id patient ID
#' @param pep_id peptide ID
#' @param start_time duration of infection of the earliest sample (in years)
#' @param end_time duration of infection of the latest sample (in years)
#'
#' @return data frame with `pt_id`, `pep_id`, `slope`, `pval`, `pep_corr`
get_rate_and_corr <- function(pt_id, pep_id, start_time, end_time){
  
  print(paste(pt_id, pep_id))
  
  # Get column of peptide, years, and ptid
  rows <- which(rc_ptid_yrs$ptid == pt_id &
                  rc_ptid_yrs$yrs_post_sero >= start_time &
                  rc_ptid_yrs$yrs_post_sero <= end_time)
  
  # Subset sample anno
  data_sub <- rc_ptid_yrs[rows, c(pep_id, "yrs_post_sero")]
  
  # Rename columns
  colnames(data_sub) <- c("peptide_rc", "yrs_post_sero")
  
  # Run lm
  lm_results <- lm(peptide_rc ~ log10(yrs_post_sero), data = data_sub)
  
  # Get p-value
  pval <- pf(summary(lm_results)$fstatistic[1],
             summary(lm_results)$fstatistic[2],
             summary(lm_results)$fstatistic[3],
             lower.tail = FALSE)
  
  # Get correlation
  pep_corr <- cor(data_sub$yrs_post_sero, data_sub$peptide_rc)
  
  
  return(data.frame(ptid = pt_id,
                    pep_id = pep_id,
                    slope = unname(lm_results$coefficients[2]),
                    pval = unname(pval),
                    pep_corr = pep_corr))
}

max_duration <- ceiling(max(sample_anno$yrs_post_sero))

slopes <- apply(expand.grid(pt_anno$ptid, pep_anno$pep_id[hiv_ind]), 
                1, function(row) get_rate_and_corr(row[1], row[2], 1/6, max_duration))
