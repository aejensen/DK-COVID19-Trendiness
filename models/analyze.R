rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(DEoptim)
library(mvtnorm)
library(rstan)
library(parallel)
library(loo)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

########################################################
# Load data
########################################################
#Downloaded from
#https://covid19.ssi.dk/overvagningsdata/download-fil-med-overvaagningdata
d <- read.csv("../data/Newly_admitted_over_time.csv", sep=";")
tObs <- as.numeric(as.Date(d$Dato) - as.Date("2020-03-01"))

dat <- list(N = nrow(d),
            t = tObs - mean(tObs), #center time - important!
            y = d$Total,
            P = 200)
dat$L <- 5/2 * max(dat$t) #boundary correction factor

########################################################
# Fit model
########################################################
iter <- 10 * 10^3

m <- stan_model("../models/gpTrend_negbinom_Matern.stan")
samps <- sampling(m, 
                  data = dat, 
                  iter = iter,
                  seed = 12345, 
                  chains = 4,
                  cores = 4)

#Posterior of hyper-parameters
hist(extract(samps, "m")$m)
hist(extract(samps, "alpha")$alpha)
hist(extract(samps, "rho")$rho)
hist(extract(samps, "nu")$nu)

mu_post <- extract(samps, "mu")$mu
dmu_post <- extract(samps, "dmu")$dmu
y_pred <- extract(samps, "yPred")$yPred
TDI <- apply(dmu_post, 2, function(q) mean(q > 0))

########################################################################
# Plots
########################################################################
band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

SAVE <- FALSE
tMax <- 224 + 14*12
yfmax <- 250
plotLandMarks <- FALSE

if(SAVE) {
  pdf("../figures/negbinom_fig1_2020.pdf", width = 16, height = 6.5)
}

par(mfrow=c(1,2), bty="n", mar = c(3, 3, 1, 0), mgp=c(2,1,0))

plot(tObs, dat$y, pch = 19, xlab="Antal dage siden 1. marts 2020", 
     ylab="Antal", type="n", ylim=c(0, 250), xaxt="n", xlim=c(0, tMax), cex.axis=0.75)
axis(1, seq(0, tMax, 14), cex.axis=0.7)
lines(rep(as.Date("2020-03-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-04-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-05-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-06-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-07-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-08-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-09-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-10-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-11-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2020-12-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2021-01-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2021-02-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")
lines(rep(as.Date("2021-03-01") - as.Date("2020-03-01"), 2), c(0, 250), lty=3, col="gray50")

band(tObs, 
     apply(y_pred, 2, quantile, 0.025),
     apply(y_pred, 2, quantile, 0.975),
     col = "gray90")
band(tObs, 
     apply(mu_post, 2, quantile, 0.025),
     apply(mu_post, 2, quantile, 0.975),
     col = "gray65")
lines(tObs, apply(mu_post, 2, mean), lwd=2)
points(tObs, dat$y, pch=19, cex=0.5)

legend(0, 250,
       c("Gennemsnit", "95% sandsynlighedsinterval", "95% prædiktionsinterval"), 
       col = c("black", "gray65", "gray90"), lwd = 2, bty="n", cex=0.7, 
       lty = c(1, NA, NA), pch = c(NA, 15, 15), pt.cex=1.5)
title("Antal daglige indlæggelser", font.main=1)

text(as.Date("2020-03-15") - as.Date("2020-03-01"), yfmax, "Mar", pos=3, cex=0.8)
text(as.Date("2020-04-15") - as.Date("2020-03-01"), yfmax, "Apr", pos=3, cex=0.8)
text(as.Date("2020-05-15") - as.Date("2020-03-01"), yfmax, "Maj", pos=3, cex=0.8)
text(as.Date("2020-06-15") - as.Date("2020-03-01"), yfmax, "Jun", pos=3, cex=0.8)
text(as.Date("2020-07-15") - as.Date("2020-03-01"), yfmax, "Jul", pos=3, cex=0.8)
text(as.Date("2020-08-15") - as.Date("2020-03-01"), yfmax, "Aug", pos=3, cex=0.8)
text(as.Date("2020-09-15") - as.Date("2020-03-01"), yfmax, "Sep", pos=3, cex=0.8)
text(as.Date("2020-10-15") - as.Date("2020-03-01"), yfmax, "Okt", pos=3, cex=0.8)
text(as.Date("2020-11-15") - as.Date("2020-03-01"), yfmax, "Nov", pos=3, cex=0.8)
text(as.Date("2020-12-15") - as.Date("2020-03-01"), yfmax, "Dec", pos=3, cex=0.8)
text(as.Date("2021-01-15") - as.Date("2020-03-01"), yfmax, "Jan", pos=3, cex=0.8)
text(as.Date("2021-02-15") - as.Date("2020-03-01"), yfmax, "Feb", pos=3, cex=0.8)
text(as.Date("2021-03-15") - as.Date("2020-03-01"), yfmax, "Mar", pos=3, cex=0.8)

plot(tObs, apply(dmu_post, 2, mean), lwd = 2, type="n", yaxt="n", xlim=c(0,tMax),
     ylim=c(-10, 10), xlab="Antal dage siden 1. marts 2020", 
     ylab="Hældning", xaxt="n")
axis(2, seq(-10, 10, 2), cex.axis=0.7)
axis(1, seq(0, tMax, 14), cex.axis=0.7)

lines(rep(as.Date("2020-03-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-04-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-05-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-06-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-07-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-08-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-09-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-10-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-11-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2020-12-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2021-01-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2021-02-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")
lines(rep(as.Date("2021-03-01") - as.Date("2020-03-01"), 2), c(-10, 10), lty=3, col="gray50")

band(tObs, 
     apply(dmu_post, 2, quantile, prob = 0.025), 
     apply(dmu_post, 2, quantile, prob = 0.975), col = "gray65")
lines(tObs, apply(dmu_post, 2, mean), lwd = 2)
abline(h = 0, lty = 2)
title("Trend for antal indlæggelser", font.main=1)

text(as.Date("2020-03-15") - as.Date("2020-03-01"), 10, "Mar", pos=3, cex=0.8)
text(as.Date("2020-04-15") - as.Date("2020-03-01"), 10, "Apr", pos=3, cex=0.8)
text(as.Date("2020-05-15") - as.Date("2020-03-01"), 10, "Maj", pos=3, cex=0.8)
text(as.Date("2020-06-15") - as.Date("2020-03-01"), 10, "Jun", pos=3, cex=0.8)
text(as.Date("2020-07-15") - as.Date("2020-03-01"), 10, "Jul", pos=3, cex=0.8)
text(as.Date("2020-08-15") - as.Date("2020-03-01"), 10, "Aug", pos=3, cex=0.8)
text(as.Date("2020-09-15") - as.Date("2020-03-01"), 10, "Sep", pos=3, cex=0.8)
text(as.Date("2020-10-15") - as.Date("2020-03-01"), 10, "Okt", pos=3, cex=0.8)
text(as.Date("2020-11-15") - as.Date("2020-03-01"), 10, "Nov", pos=3, cex=0.8)
text(as.Date("2020-12-15") - as.Date("2020-03-01"), 10, "Dec", pos=3, cex=0.8)
text(as.Date("2021-01-15") - as.Date("2020-03-01"), 10, "Jan", pos=3, cex=0.8)
text(as.Date("2021-02-15") - as.Date("2020-03-01"), 10, "Feb", pos=3, cex=0.8)
text(as.Date("2021-03-15") - as.Date("2020-03-01"), 10, "Mar", pos=3, cex=0.8)

if(SAVE) {
  dev.off()
}

########################################################################
# Figure 2
########################################################################
if(SAVE) {
  pdf("../figures/negbinom_fig2_2020.pdf", width = 10, height = 5)
}

par(mfrow=c(1,1), bty="n", mar = c(2.5, 2.3, 1, 0), mgp=c(1.4,0.4,0))
plot(tObs, TDI*100, type="n", lty = 1, lwd = 2,
     xlab="Antal dage siden 1. marts 2020", 
     ylab="Sandsynlighed for voksende antal indlæggelser [%]", 
     ylim=c(0,100), xaxt="n", xlim=c(0, tMax))

if(plotLandMarks) {
  lines(rep(as.Date("2020-03-11") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
  points(as.Date("2020-03-11") - as.Date("2020-03-01"), 80, pch=21, bg="white")
  lines(rep(as.Date("2020-04-15") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
  points(as.Date("2020-04-15") - as.Date("2020-03-01"), 80, pch=21, bg="white")
  lines(rep(as.Date("2020-04-20") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
  points(as.Date("2020-04-20") - as.Date("2020-03-01"), 80, pch=21, bg="white")
  lines(rep(as.Date("2020-05-18") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
  points(as.Date("2020-05-18") - as.Date("2020-03-01"), 80, pch=21, bg="white")
  lines(rep(as.Date("2020-05-27") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
  points(as.Date("2020-05-27") - as.Date("2020-03-01"), 80, pch=21, bg="white")

  text(as.Date("2020-03-11") - as.Date("2020-03-01"), 80, "Danmark\nlukker", pos = 3, cex=0.7)
  text(as.Date("2020-04-15") - as.Date("2020-03-01"), 80, "Skoler op til\n5. klasse åbner", pos = 3, cex=0.7)
  text(as.Date("2020-04-20") - as.Date("2020-03-01"), 87, "Frisører og køre-\nskoler åbner", pos = 3, cex=0.7)
  text(as.Date("2020-05-18") - as.Date("2020-03-01"), 80, "Restauranter\nåbner", pos = 3, cex=0.7)
  text(as.Date("2020-05-27") - as.Date("2020-03-01"), 83, "Gymnasier\nåbner", pos = 3, cex=0.7)
}

lines(rep(as.Date("2020-04-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-05-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-06-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-07-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-08-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-09-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-10-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-11-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-12-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2021-01-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")

lines(tObs, TDI*100, type="l", lty = 1, lwd = 2)
axis(1, seq(0, tMax, 14), cex.axis=0.90)
abline(h = 50, lty = 2)
title("Trend Direction Index", font.main=1)

text(as.Date("2020-03-15") - as.Date("2020-03-01"), 99, "Marts", pos=3, cex=0.8)
text(as.Date("2020-04-15") - as.Date("2020-03-01"), 99, "April", pos=3, cex=0.8)
text(as.Date("2020-05-15") - as.Date("2020-03-01"), 99, "Maj", pos=3, cex=0.8)
text(as.Date("2020-06-15") - as.Date("2020-03-01"), 99, "Juni", pos=3, cex=0.8)
text(as.Date("2020-07-15") - as.Date("2020-03-01"), 99, "Juli", pos=3, cex=0.8)
text(as.Date("2020-08-15") - as.Date("2020-03-01"), 99, "August", pos=3, cex=0.8)
text(as.Date("2020-09-15") - as.Date("2020-03-01"), 99, "September", pos=3, cex=0.8)
text(as.Date("2020-10-15") - as.Date("2020-03-01"), 99, "Oktober", pos=3, cex=0.8)
text(as.Date("2020-11-15") - as.Date("2020-03-01"), 99, "November", pos=3, cex=0.8)
text(as.Date("2020-12-15") - as.Date("2020-03-01"), 99, "December", pos=3, cex=0.8)

points(max(tObs), TDI[length(tObs)]*100, pch=19, cex=0.6)

text(max(tObs), round(TDI[length(tObs)]*100, 2), 
     paste(round(TDI[length(tObs)]*100, 2), "%", sep=""), pos=1, cex=0.7)

if(SAVE) {
	dev.off()
}

