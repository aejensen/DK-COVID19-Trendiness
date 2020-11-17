/* 
  Latent Gaussian Process regression model with negative binomial outcome model
  using the Hilbert space reduced rank method of Solin and Särkkä (2020) and 
  Riutort-Mayol, Bürkner, Andersen, Solin and Vehtari (2020). 
   
  Andreas Kryger Jensen, 2020. 
*/

#include /gpTrend_approx_functions.stan

data {
  int<lower=1> N;   //number of observation points
  vector[N] t;      //time values
  int y[N];         //integer outcome
  
  real L;           //boundary condition factor
  int<lower=1> P;   //number of basis functions
}

transformed data {
  /* Pre-calculate matrices of eigenvectors and derivatives */
  matrix[N, P] PHI;
  matrix[N, P] PHI_d1;
  matrix[N, P] PHI_d2;
  
  for (p in 1:P) { 
    PHI[, p] = phi(L, p, t); 
    PHI_d1[, p] = phi_d1(L, p, t);
    PHI_d2[, p] = phi_d2(L, p, t);
  }
}

parameters {
	vector[P] beta;        //coefficients for f

  real m;                //mean parameter
  real<lower = 0> alpha; //variance parameter
  real<lower = 0> rho;   //length scale parameter

  real<lower = 0> theta; //dispersion parameter
}

transformed parameters{
  vector[N] f;
  vector[P] diagSPD;  //"diagonal" Delta matrix with entries S_theta(lambda^1/2)
  vector[P] SPD_beta; //Delta_PxP * beta_1xP

  for(p in 1:P){ 
    diagSPD[p] = sqrt(spd_SE(alpha, rho, sqrt(lambda(L, p)))); 
  }
   
  //S(sqrt(lambda_p))^(1/2) * beta_p, p = 1,..., P
  SPD_beta = diagSPD .* beta;
	  
  /* f = m + PHI_(n x P) * [S(sqrt(lambda_1))^(1/2) * beta_1, ...]_(P x 1)
       = m + \sum_{p=1}^P Phi_p * S(sqrt(lambda_p))^(1/2) * beta_p */
  f = m + PHI * (diagSPD .* beta);
}

model{
  m ~ normal(0, 3);

  alpha ~ normal(0, 3);
  rho ~ normal(0, 5);
  beta ~ normal(0, 1);

  theta ~ gamma(0.01, 0.01);

  y ~ neg_binomial_2_log(f, theta);
}

generated quantities {
  vector[N] yPred;    //posterior predictive values
  vector[N] log_lik;  //log likelihood values
  
  vector[N] df = PHI_d1 * SPD_beta;
  vector[N] ddf = PHI_d2 * SPD_beta;
  
  /* Calculate posterior mean and first and second derivatives using the chain rule */
  vector[N] mu = exp(f);
  vector[N] dmu = mu .* df;
  vector[N] ddmu = dmu .* df - mu .* ddf;
  
  for(i in 1:N) {
    /* truncate due to possible numerical overflow */
    yPred[i] = f[i] > 20 ? neg_binomial_2_log_rng(20, theta)
                         : neg_binomial_2_log_rng(f[i], theta);
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | f[i], theta);
  }
}
