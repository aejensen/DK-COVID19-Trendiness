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
#https://www.ssi.dk/sygdomme-beredskab-og-forskning/sygdomsovervaagning/c/covid19-overvaagning
d <- read.csv("../data/Newly_admitted_over_time.csv", sep=";")
tObs <- as.numeric(as.Date(d$Dato) - as.Date("2020-03-01"))

dat <- list(N = nrow(d),
            t = tObs - mean(tObs), #center time - important!
            y = d$Total,
            L = 5/2 * max(tObs),
            P = 125)

########################################################
# Fit models
########################################################
iter <- 1 * 10^4

m_negbinom_SE <- stan_model("../models/gpTrend_negbinom_SE.stan")
samp_negbinom_SE <- sampling(m_negbinom_SE, 
                             data = dat, 
                             iter = iter, 
                             seed = 12345, 
                             chains = 4, 
                             cores = 4)

m_negbinom_Matern <- stan_model("../models/gpTrend_negbinom_Matern.stan")
samp_negbinom_Matern <- sampling(m_negbinom_Matern, 
                                 data = dat, 
                                 iter = iter,
                                 seed = 12345, 
                                 chains = 4,
                                 cores = 4)

m_normal_SE <- stan_model("../models/gpTrend_normal_SE.stan")
samp_norm_SE <- sampling(m_normal_SE, 
                         data = dat, 
                         iter = iter, 
                         seed = 12345,
                         chains = 4, 
                         cores = 4)


########################################################
# Compare models
########################################################
getLOO <- function(fit) {
  log_lik_1 <- extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- relative_eff(exp(log_lik_1), cores = 4) 
  loo(log_lik_1, r_eff = r_eff, cores = 4)
}

loo_negbinom_SE <- getLOO(samp_negbinom_SE)
loo_negbinom_Matern <- getLOO(samp_negbinom_Matern)
loo_norm_SE <- getLOO(samp_norm_SE)

loo_compare(loo_negbinom_SE, loo_negbinom_Matern, loo_norm_SE)
waic(loo_negbinom_SE)

plot(loo_negbinom_Matern)
plot(loo_negbinom_SE)
plot(loo_norm_SE)

########################################################
# Posterior diagnostics
########################################################
hist(extract(samp_negbinom_Matern, "nu")$nu)


