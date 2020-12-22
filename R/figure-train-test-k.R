# figure-train-test.R ----------
# Code to generate:
#   * `figures/train-test-k.pdf`
# -----------------------
source(here("R/run-ktsp.R"))

train_accuracies <- train_ktsp %>%
  group_by(num_pairs) %>%
  summarize(accuracy = mean(recent == kpred), 
            .groups = "drop")

test_accuracies <- test_ktsp %>%
  group_by(num_pairs) %>%
  summarize(accracy = mean(recent == kpred), 
            .groups = "drop")

pdf("figures/train-test-k.pdf",width=7,height=7)
par(las=1,cex.axis=1.2,cex.lab=1.2,yaxs="i")
plot(c(1,10),c(0.8,1),type="n",xaxt="n",yaxt="n",
     xlab="k (number of peptides in the classifier)",ylab="accuracy")
abline(v=1:10,lty="dotted",col="lightgrey")
abline(h=seq(0.8,1,0.025),lty="dotted",col="lightgrey")
lines(train_accuracies$accuracy,type="b",col="blue",cex=1.3,lwd=2)
lines(test_accuracies$accracy,type="b",col="red",cex=1.3,lwd=2)
axis(1,1:10)
axis(2,seq(0.80,1,0.05))
legend("bottomright",c("training","test"),lwd=2,
       col=c("blue","red"),bg="white",cex=1.3)
dev.off()



