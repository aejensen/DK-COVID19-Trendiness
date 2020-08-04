rm(list=ls())

library(rstan)
library(parallel)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m <- stan_model("gpNegBinomApprox.stan")

# Fix a bug in parallel, RStudio and R > 4.0.0 on Mac
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

########################################################
# Load data
########################################################
d <- read.csv("../data/Newly_admitted_over_time.csv", sep=";")
tObs <- as.numeric(as.Date(d$Dato) - as.Date("2020-03-01"))
#tPred <- seq(min(tObs), max(tObs), length.out = 2*length(tObs))

########################################################
# Fit model
########################################################
sDat <- list(N = nrow(d),
             t = tObs - mean(tObs), #center time - important!
             y = d$Total,
             L = 5/2 * max(tObs),
             P = 120)

samps <- sampling(m, data = sDat, iter = 10000, seed = 12345, chains = 4, cores = 4)

save.image(file = "DK-COVID19_negbinom.RData") #too big for GitHub

#Trace plot for hyper-parameters
#mcmc_trace(samps,  pars = c("m", "alpha", "rho", "theta"), n_warmup = 25000/2)

########################################################
#
########################################################
load("DK-COVID19_negbinom.RData")

mu_post <- extract(samps, "mu")$mu
y_pred <- extract(samps, "yPred")$yPred

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("../figures/negbinom_fig1.pdf", width = 8, height = 5)
par(mfrow=c(1,1), bty = "n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))

plot(tObs, sDat$y, ylim=c(0,120), type="n",
     xaxt="n", xlab="Number of days since March 1th 2020",
     ylab="Number", xlim=c(0,168))
axis(1, seq(0, 168, 14), cex.axis=0.95)
band(tObs, 
     apply(y_pred, 2, quantile, 0.025),
     apply(y_pred, 2, quantile, 0.975),
     col = "gray90")
band(tObs, 
     apply(mu_post, 2, quantile, 0.025),
     apply(mu_post, 2, quantile, 0.975),
     col = "gray65")
lines(tObs, apply(mu_post, 2, mean), lwd=2)
points(tObs, sDat$y, pch=19, cex=0.5)

title("Daily hospital admissions", font.main=1)

legend("topleft",
       c("Posterior mean", "95% CI", "95% posterior prediction interval"), 
       col = c("black", "gray65", "gray90"), lwd = 2, bty="n", cex=0.7, 
       lty = c(1, NA, NA), pch = c(NA, 15, 15), pt.cex=1.5)
dev.off()

df <- t(apply(mu_post, 1, function(f) {
  predict(smooth.spline(tObs, f, all.knots = TRUE), deriv=1)$y
}))

tdi <- apply(df, 2, function(q) mean(q > 0))


pdf("../figures/negbinom_fig2.pdf", width = 12, height = 6)
par(mfrow=c(1, 2), bty = "n", mar = c(2.3, 2.3, 1, 0.5), mgp=c(1.3,0.4,0))

plot(tObs, sDat$y, ylim=c(-6,8), type="n",
     xaxt="n", xlab="Antal dage siden 1. marts 2020",
     ylab="Hælding")
axis(1, seq(0, 168, 14), cex.axis=0.90)
band(tObs, 
     apply(df, 2, quantile, 0.025),
     apply(df, 2, quantile, 0.975),
     col = "gray65")
lines(tObs, apply(df, 2, mean), lwd=2)
abline(h = 0, lty=2)

plot(tObs, tdi * 100, ylim=c(0,100), type="l",
     xaxt="n", xlab="Antal dage siden 1. marts 2020",
     ylab="Trend Direction Index [%]", lwd=2, xlim=c(0,168))
axis(1, seq(0, 168, 14), cex.axis=0.90)
abline(h = 50, lty=2)
dev.off()

