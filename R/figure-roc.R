# figure-roc.R ----------
# Code to generate:
#   * `figures/roc.pdf`
#   * `figures/prc.pdf`
#   * `figures/roc_with_window.pdf`
# -----------------------

library(RColorBrewer)
library(tidyverse)

source(here("R/run-ktsp.R"))

train_opti_sens <- train_curve %>% 
  filter(pep_pair %in% c("pep_25409-pep_77241",
                         "pep_93864-pep_92687",
                         "pep_26777-pep_20904",
                         "pep_77323-pep_32119")) %>%
  left_join(train_opti_cutoffs %>% 
              dplyr::select(pep_pair, opti_cutoff), by = "pep_pair") %>%
  filter(cutoff == opti_cutoff)

pps <- c("pep_25409-pep_77241","pep_93864-pep_92687",
         "pep_26777-pep_20904","pep_77323-pep_32119")

ppsn <- pps
ppsn <- gsub("pep_","",ppsn)
ppsn <- gsub("-"," / ",ppsn)

dat <- subset(train_curve, pep_pair%in%pps)
dat$pep_pair <- factor(dat$pep_pair)
datlag <- lag_curve[[1]]

mcol <- brewer.pal(9,"Blues")[c(4,6,8,9)]
mcol <- brewer.pal(10,"RdBu")[7:10]

# PRC -----------------------
pdf("figures/prc.pdf",width=7,height=7)
par(pty="s",bty="n",xaxs="i",yaxs="i",las=1,mgp=c(3.5,1,0),
    font.lab=2,cex.lab=1.2)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",
     xlab="Recall  ( Sensitivity )",
     ylab="Precision  ( Positive Predictive Value )")
abline(v=c(0,1))
abline(h=c(0,1))
abline(h=seq(0,1,0.1),lty="dotted",col="lightgrey")
abline(v=seq(0,1,0.1),lty="dotted",col="lightgrey")
lines(datlag$sens,datlag$ppv,col="red",lwd=2)
for(i in 1:4){
  tmp <- subset(dat,pep_pair==pps[i])
  lines(tmp$sens,tmp$ppv,col=mcol[i],lwd=1.5)
}
wh <- c(1,max(which(datlag$ppv==1)))
lines(datlag$sens[wh],rep(1,2),col="red",lwd=3)
wh <- c(min(which(datlag$sens==1)),nrow(datlag))
lines(rep(1,2),datlag$ppv[wh],col="red",lwd=4)
for(i in 1:4){
  tmp <- subset(train_opti_sens,pep_pair==pps[i])
  points(tmp$sens,tmp$ppv,pch=20,col=mcol[i],cex=2)
}
axis(1,seq(0,1,0.1),font.axis=2)
axis(2,seq(0,1,0.1),font.axis=2)
nms <- c(paste0("Pair ",1:4,"   "),"LAg")
legend("bottomright",nms,lty=1,lwd=2,col=c(mcol,"red"),bg="white",cex=1.0)
dev.off()

# ROC -----------------------
pdf("figures/roc.pdf",width=7,height=7)
par(pty="s",bty="n",xaxs="i",yaxs="i",las=1,mgp=c(3.5,1,0),
    font.lab=2,cex.lab=1.2)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",
     xlab="False Positive Rate  ( 1 - Specificity )",
     ylab="True Positive Rate  ( Sensitivity )")
abline(v=c(0,1))
abline(h=c(0,1))
abline(h=seq(0,1,0.1),lty="dotted",col="lightgrey")
abline(v=seq(0,1,0.1),lty="dotted",col="lightgrey")
lines(1-datlag$spec,datlag$sens,col="red",lwd=2)
auc <- NULL
for(i in 1:4){
  tmp <- subset(dat,pep_pair==pps[i])
  lines(1-tmp$spec,tmp$sens,col=mcol[i],lwd=1.5)
  auc <- c(auc,tmp$auc_roc[1])
}
for(i in 1:4){
  tmp <- subset(train_opti_sens,pep_pair==pps[i])
  points(1-tmp$spec,tmp$sens,pch=20,col=mcol[i],cex=2)
}
auc <- c(auc,datlag$auc_roc[1])
axis(1,seq(0,1,0.1),font.axis=2)
axis(2,seq(0,1,0.1),font.axis=2)
auc <- paste0("( ",formatC(auc,digits=2,format="f")," )")
nms <- c(paste0("Pair ",1:4,"   "),"LAg")
legend("bottomright",paste(auc,nms,sep="  "),
       lty=1,lwd=2,col=c(mcol,"red"),bg="white",cex=1.0)
dev.off()

# ROC with window -----------------------
pep_train_window <- pep_ratios %>%
  filter(ptid %in% ptid_train &
           yrs_post_sero >= 1/6) %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0))

training_window_cutoffs <- lapply(rownames(pairs), function(x) get_cutoffs(pep_train_window, x)) 

# data for ROC/PRC curves
train_window_curve <- training_window_cutoffs %>%
  lapply(function(x) x[[1]]) %>%
  plyr::ldply()

# df of optimal cutofs for each peptide pairs
train_window_opti_cutoffs <- training_window_cutoffs %>%
  lapply(function(x) x[[2]]) %>%
  plyr::ldply()

train_window_opti_sens <- train_window_curve %>% 
  filter(pep_pair %in% c("pep_25409-pep_77241",
                         "pep_93864-pep_92687",
                         "pep_26777-pep_20904",
                         "pep_77323-pep_32119")) %>%
  left_join(train_window_opti_cutoffs %>% 
              dplyr::select(pep_pair, opti_cutoff), by = "pep_pair") %>%
  filter(cutoff == opti_cutoff)

# lag data
lag_window_train <- sample_anno %>%
  filter(yrs_post_sero >= 1/6 &
           (ptid %in% ptid_train)) %>%
  mutate(recent = ifelse(yrs_post_sero <= 1, 1, 0)) %>%
  dplyr::select(ptid, yrs_post_sero, recent, lag)

lag_window_curve <- get_cutoffs(lag_window_train, "lag")

dat <- subset(train_window_curve, pep_pair%in%pps)
dat$pep_pair <- factor(dat$pep_pair)
datlag <- lag_window_curve[[1]]

pdf("figures/roc_with_window.pdf",width=7,height=7)
par(pty="s",bty="n",xaxs="i",yaxs="i",las=1,mgp=c(3.5,1,0),
    font.lab=2,cex.lab=1.2)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",
     xlab="False Positive Rate  ( 1 - Specificity )",
     ylab="True Positive Rate  ( Sensitivity )")
abline(v=c(0,1))
abline(h=c(0,1))
abline(h=seq(0,1,0.1),lty="dotted",col="lightgrey")
abline(v=seq(0,1,0.1),lty="dotted",col="lightgrey")
lines(1-datlag$spec,datlag$sens,col="red",lwd=2)
auc <- NULL
for(i in 1:4){
  tmp <- subset(dat,pep_pair==pps[i])
  lines(1-tmp$spec,tmp$sens,col=mcol[i],lwd=1.5)
  auc <- c(auc,tmp$auc_roc[1])
}
for(i in 1:4){
  tmp <- subset(train_window_opti_sens,pep_pair==pps[i])
  points(1-tmp$spec,tmp$sens,pch=20,col=mcol[i],cex=2)
}
auc <- c(auc,datlag$auc_roc[1])
axis(1,seq(0,1,0.1),font.axis=2)
axis(2,seq(0,1,0.1),font.axis=2)
auc <- paste0("( ",formatC(auc,digits=2,format="f")," )")
nms <- c(paste0("Pair ",1:4,"   "),"LAg")
legend("bottomright",paste(auc,nms,sep="  "),
       lty=1,lwd=2,col=c(mcol,"red"),bg="white",cex=1.0)
dev.off()
