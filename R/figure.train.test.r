library(tidyverse)
load("data/tsp.rda")

# Get accuracies for all k's for test set
test_tsp_votes <- pep_ratio_sm_post2mo_test %>%
  gather("pep_pair", "rc", -ptid, -yrs_post_sero, -recent) %>%
  left_join(train_window_summary, by = "pep_pair") %>%
  mutate(pred_rec = ifelse(rc<= opti_cutoff, 1, 0)) %>%
  dplyr::select(ptid, yrs_post_sero,  pep_pair, recent, pred_rec) %>%
  filter(pep_pair %in% kpep_list$pep_pair) %>%
  spread(pep_pair, pred_rec)

# Count votes
k_test_votes <- apply(test_tsp_votes[,kpep_list$pep_pair], 1, cumsum) %>% t()
rownames(k_test_votes) <- paste0(test_tsp_votes$ptid, "_", test_tsp_votes$yrs_post_sero)
colnames(k_test_votes) <- paste0("k_", 1:10)

ktsp_test_outcomes <- k_test_votes/matrix(rep(1:10, nrow(k_test_votes)), ncol = 10, byrow = TRUE)
ktsp_test_outcomes <- ifelse(ktsp_test_outcomes > 0.5, 1, 0)

ktsp_test_accurate <- ktsp_test_outcomes - matrix(rep(test_tsp_votes$recent, 10), ncol = 10)
ktsp_test_accurate <- ifelse(ktsp_test_accurate == 0, 1, 0) 

k_test_window <- data.frame(ktsp_test_accurate) %>%
  mutate(sample = rownames(ktsp_test_accurate)) %>%
  separate(sample, c("ptid", "yrs_post_sero"), sep = "_") %>%
  filter(yrs_post_sero >= 1.5 | yrs_post_sero <= 0.5)

k_test_wind_accuracy <- apply(as.matrix(k_test_window[, paste0("k_", 1:10)]), 2, mean)

pdf("plots/traintest.pdf",width=7,height=7)
par(las=1,cex.axis=1.2,cex.lab=1.2,yaxs="i")
plot(c(1,10),c(0.8,1),type="n",xaxt="n",yaxt="n",
     xlab="k (number of peptides in the classifier)",ylab="accuracy")
abline(v=1:10,lty="dotted",col="lightgrey")
abline(h=seq(0.8,1,0.025),lty="dotted",col="lightgrey")
lines(k_wind_accuracy,type="b",col="blue",cex=1.3,lwd=2)
lines(k_test_wind_accuracy,type="b",col="red",cex=1.3,lwd=2)
axis(1,1:10)
axis(2,seq(0.80,1,0.05))
legend("bottomright",c("training","test"),lwd=2,
       col=c("blue","red"),bg="white",cex=1.3)
dev.off()



