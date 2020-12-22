# figure-loocv-k.R ----------
# Code to generate:
#   * `figures/loocv-k.pdf`
# -----------------------

source(here("R/run-cv.R"))

dat <- matrix(cv_train$accuracy,ncol=10,byrow=TRUE)[,c(1,3:10,2)]
cvtr <- apply(dat,2,mean)

dat <- matrix(cv_test$accuracy,ncol=10,byrow=TRUE)[,c(1,3:10,2)]
cvte <- apply(dat,2,mean)

pdf("figures/loocv-k.pdf",width=7,height=7)
par(las=1,cex.axis=1.2,cex.lab=1.2,yaxs="i")
plot(c(1,10),c(0.8,1),type="n",xaxt="n",yaxt="n",
     xlab="k (number of peptides in the classifier)",ylab="accuracy")
abline(v=1:10,lty="dotted",col="lightgrey")
abline(h=seq(0.8,1,0.025),lty="dotted",col="lightgrey")
lines(cvtr,type="b",col="blue",cex=1.3,lwd=2)
lines(cvte,type="b",col="red",cex=1.3,lwd=2)
axis(1,1:10)
axis(2,seq(0.80,1,0.05))
legend("bottomright",c("leave-one-out training","leave-one-out evaluation"),lwd=2,
       col=c("blue","red"),bg="white",cex=1.3)
dev.off()
