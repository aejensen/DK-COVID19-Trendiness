rm(list=ls())

library(DEoptim)
library(mvtnorm)
library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Fix a bug in parallel, RStudio and R > 4.0.0 on Mac
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

########################################################
# Load data
########################################################
#Downloaded from
#https://www.ssi.dk/sygdomme-beredskab-og-forskning/sygdomsovervaagning/c/covid19-overvaagning
d <- read.csv("../data/Newly_admitted_over_time.csv", sep=";")

tObs <- as.numeric(as.Date(d$Dato) - as.Date("2020-03-01"))
tPred <- seq(min(tObs), max(tObs), length.out = ceiling(length(tObs)*1))

########################################################
# Fit model
########################################################
optPar <- function(t, y) {
  rqCov <- function(s, t, alpha, rho, nu) {
    alpha^2 * (1 + (s-t)^2 / (2 * nu * rho^2))^(-nu)
  }
  
  ctl <- DEoptim.control(itermax = 1000, trace = 100)
  set.seed(1234)
  opt <- DEoptim(function(par) {
    mu <- rep(par[1], length(t))
    cMat <- outer(t, t, rqCov, par[2], par[3], par[4]) + diag(par[5]^2 , length(t))
    -mvtnorm::dmvnorm(y, mu, cMat, log=TRUE)
  }, lower = c(0,0,0,0,0), upper = c(max(y),100,50,1000,50), control = ctl)
  opt$optim$bestmem
}

pars <- optPar(tObs, d$Total)

sDat <- list(n = length(tObs), t = tObs, y = d$Total, p = length(tPred), tPred = tPred)
sDat$mu <- pars[1]
sDat$alpha <- pars[2]
sDat$rho <- pars[3]
sDat$nu <- pars[4]
sDat$sigma <- pars[5]

m <- stan_model("gptrendFixed.stan")

fit <- sampling(m, data = sDat, iter = 10000, seed = 12345, algorithm = "Fixed_param", cores = 4)
posterior <- extract(fit, "pred")$pred

save.image("DK-COVID19.RData") #too big for GitHub
