# helper-functions.R

# Create peptide ratios -----------------------
#' Function to calculate log10 of relative ratios of the normalized RC for each peptide pair
#'
#' @param pep1 upgoing peptide
#' @param pep2 downgoing peptide
#' @param data data set that has `ptid`, `yrs_post_sero`, and the peptides in the columns
#'
#' @return datafram with `ptid`, `yrs_post_sero`, `pep_pair`, `logratio`
get_ratio <- function(pep1, pep2, data){
  
  logratio <- unname(data[, pep1] - data[, pep2])
  
  return(data.frame(ptid = data$ptid,
                    yrs_post_sero = data$yrs_post_sero,
                    pep_pair = paste(pep1, pep2, sep = "-"),
                    logratio = logratio))
}
