# figure-logistic.R ----------
# Code to generate:
#   * `figures/logistic.pdf`
# -----------------------

library(RColorBrewer)
library(tidyverse)

source(here("R/run-ktsp.R"))

# Get lag predictions for all samples after 2 months
lag_pred <- sample_anno %>%
  filter(yrs_post_sero >= 1/6) %>%
  mutate(maa_rec = ifelse(lag >= 1.5 | (lag < 1.5 & viral_load <= 3), 0, 1)) %>%
  select(key, ptid, yrs_post_sero, lag, viral_load, maa_rec)

# Get current kTSP classification for all samples after 2 months
ktsp_pred <- pep_ratios %>%
  select(ptid, yrs_post_sero, kpep_list$pep_pair[1:4]) %>%
  pivot_longer(cols = contains("pep"), 
               names_to = "pep_pair", 
               values_to = "rc") %>%
  left_join(train_opti_cutoffs[, c("pep_pair", "opti_cutoff")], by = "pep_pair") %>%
  mutate(pred_rec = ifelse(rc <= opti_cutoff, 1, 0)) %>%
  select(ptid, yrs_post_sero, pep_pair, pred_rec) %>%
  group_by(ptid, yrs_post_sero) %>%
  summarize(rec_votes = sum(pred_rec), .groups = "drop") %>%
  mutate(ktsp_rec = ifelse(rec_votes >= 3, 1, 0)) 

# Merge the two data frames 
comparison_pred <- lag_pred %>%
  left_join(ktsp_pred[, c("ptid", "yrs_post_sero", "ktsp_rec")], by = c("ptid", "yrs_post_sero")) %>%
  mutate(acc_rec = ifelse(yrs_post_sero <= 1, 1, 0)) 

# Logistic regression fits for current 
lag_logistic_fit <- glm(maa_rec ~ yrs_post_sero,
                       data = comparison_pred,
                       family = binomial(link = "logit"))

kstp_logistic_fit <- glm(ktsp_rec ~ yrs_post_sero,
                        data = comparison_pred,
                        family = binomial(link = "logit"))

hm <- 101
xx <- seq(0,2,length=hm)
prmaa <- predict(lag_logistic_fit, data.frame(yrs_post_sero=xx))
prmaa <- exp(prmaa)/(1+exp(prmaa))
prktsp <- predict(kstp_logistic_fit, data.frame(yrs_post_sero=xx))
prktsp <- exp(prktsp)/(1+exp(prktsp))

w <- xx[2]-xx[1]
auc.maa <- sum((prmaa[-1]+prmaa[-hm])/2*w)
auc.ktsp <- sum((prktsp[-1]+prktsp[-hm])/2*w)
auc <- c(auc.ktsp,auc.maa)
auc <- paste0("( ",formatC(auc,digits=2,format="f")," )")

# CIs
p.maa <- predict(lag_logistic_fit,data.frame(yrs_post_sero=xx),
                   type="link",se.fit=TRUE)
p.maa.up <- p.maa$fit+1.96*p.maa$se.fit
p.maa.lo <- p.maa$fit-1.96*p.maa$se.fit
p.maa.up <- exp(p.maa.up)/(1+exp(p.maa.up))
p.maa.lo <- exp(p.maa.lo)/(1+exp(p.maa.lo))
auc.maa.up <- sum((p.maa.up[-1]+p.maa.up[-hm])/2*w)
auc.maa.lo <- sum((p.maa.lo[-1]+p.maa.lo[-hm])/2*w)
print(round(c(auc.maa,auc.maa.lo,auc.maa.up),2))
print(round(c(auc.maa,auc.maa.lo,auc.maa.up)*365))

p.ktsp <- predict(kstp_logistic_fit,data.frame(yrs_post_sero=xx),
                  type="link",se.fit=TRUE)
p.ktsp.up <- p.ktsp$fit+1.96*p.ktsp$se.fit
p.ktsp.lo <- p.ktsp$fit-1.96*p.ktsp$se.fit
p.ktsp.up <- exp(p.ktsp.up)/(1+exp(p.ktsp.up))
p.ktsp.lo <- exp(p.ktsp.lo)/(1+exp(p.ktsp.lo))
auc.ktsp.up <- sum((p.ktsp.up[-1]+p.ktsp.up[-hm])/2*w)
auc.ktsp.lo <- sum((p.ktsp.lo[-1]+p.ktsp.lo[-hm])/2*w)
# print(round(c(auc.ktsp.lo,auc.ktsp.up),2))
# print(round(c(auc.ktsp,auc.ktsp.lo,auc.ktsp.up)*365))

# Plot
pdf("figures/logistic.pdf",width=8,height=6)
par(las=1,xaxs="i",yaxs="i",font.lab=2,cex.lab=1.2)
plot(range(xx),c(0,1),type="n",xaxt="n",yaxt="n",
     xlab="Duration of Infection  [ years ]",
     ylab="Probability of Appearing Recent")
abline(h=seq(0,1,0.1),lty="dotted",col="lightgrey")
abline(v=seq(0,10,0.25),lty="dotted",col="lightgrey")
polygon(c(xx,rev(xx)),c(p.ktsp.lo,rev(p.ktsp.up)),col=rgb(0,0,1,0.3),border=NA)
polygon(c(xx,rev(xx)),c(p.maa.lo,rev(p.maa.up)),col=rgb(1,0,0,0.3),border=NA)
lines(xx,prktsp,lwd=3,col="blue")
lines(xx,prmaa,lwd=3,col="red")
axis(1,seq(0,10,0.25),rep("",length(seq(0,10,0.25))),tcl=-0.2)
axis(1,seq(0,10,0.5),font.axis=2)
axis(2,seq(0,1,0.1),font.axis=2)
text(1.00,0.945,"4-TSP",col="blue",cex=1.3,pos=4)
text(1.00,0.845,
     paste0("mean window: ",round(auc.ktsp*365)," days"),
     col="blue",cex=1.1,pos=4)
text(1.00,0.745,
     paste0("95% confidence interval: ",round(auc.ktsp.lo*365)," - ",
            round(auc.ktsp.up*365)," days"),
     col="blue",cex=0.8,pos=4)
text(1.00,0.645,"LAg + VL",col="red",cex=1.3,pos=4)
text(1.00,0.545,
     paste0("mean window: ",round(auc.maa*365)," days"),
     col="red",cex=1.1,pos=4)
text(1.00,0.445,
     paste0("95% confidence interval: ",round(auc.maa.lo*365)," - ",
            round(auc.maa.up*365)," days"),
     col="red",cex=0.8,pos=4)
dev.off()
