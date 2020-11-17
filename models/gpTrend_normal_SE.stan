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
  vector[N] y;      //outcome
  
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

  real<lower = 0> sigma; //residual variance
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

  sigma ~ cauchy(0, 2.5);

  y ~ normal(f, sigma);
}

generated quantities {
  vector[N] yPred;    //posterior predictive values
  vector[N] log_lik;  //log likelihood values
  
  vector[N] df = PHI_d1 * SPD_beta;
  vector[N] ddf = PHI_d2 * SPD_beta;
  
  for(i in 1:N) {
    yPred[i] = normal_rng(f[i], sigma);
    log_lik[i] = normal_lpdf(y[i] | f[i], sigma);
  }
}
