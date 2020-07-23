rm(list=ls())

library(DEoptim)
library(mvtnorm)
library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m <- stan_model("gptrendFixed.stan")

# Fix a bug in parallel, RStudio and R > 4.0.0 on Mac
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

########################################################
# Load data
########################################################
#Downloaded July 23th 2020 from
#https://www.ssi.dk/sygdomme-beredskab-og-forskning/sygdomsovervaagning/c/covid19-overvaagning
d <- read.csv("../data/Newly_admitted_over_time.csv", sep=";")

tObs <- as.numeric(as.Date(d$Dato) - as.Date("2020-03-01"))
tPred <- seq(min(tObs), max(tObs), length.out = ceiling(length(tObs)*2))

########################################################
# Fit model
########################################################
getTheTrend <- function(t, y) {
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
  par <- opt$optim$bestmem
  
  sDat <- list(n = length(t), t = t, y = y, p = length(tPred), tPred = tPred)
  sDat$mu <- par[1]
  sDat$alpha <- par[2]
  sDat$rho <- par[3]
  sDat$nu <- par[4]
  sDat$sigma <- par[5]
  
  fit <- sampling(m, data = sDat, iter = 10000, seed = 12345, algorithm = "Fixed_param", cores = 4)
  posterior <- extract(fit, "pred")$pred
  
  list(t = t, y = y, par = par, post = posterior)
}

total <- getTheTrend(tObs, d$Total)

save.image("DK-COVID19.RData") #too big for GitHub
