rm(list=ls())

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m <- stan_model("gpNegBinomApprox2.stan")

.Li# Fix a bug in parallel, RStudio and R > 4.0.0 on Mac
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

########################################################
# Load data
########################################################
d <- read.csv("../data/Newly_admitted_over_time.csv", sep=";")
tObs <- as.numeric(as.Date(d$Dato) - as.Date("2020-03-01"))

########################################################
# Fit model
########################################################
sDat <- list(N = nrow(d),
             t = tObs - mean(tObs), #center time - important!
             y = d$Total,
             L = 5/2 * max(tObs),
             P = 125)

samps <- sampling(m, data = sDat, iter = 5000, seed = 12345, chains = 2, cores = 2)


install.packages("rstan")

mu_post <- extract(samps, "mu")$mu


