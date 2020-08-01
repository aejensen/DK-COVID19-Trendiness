functions {
	real lambda(real L, int p) {
		real lam = ((p*pi())/(2*L))^2;
		return lam;
	}
	
	real spd(real alpha, real rho, real w) {
		real S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
		return S;
	}
	
	vector phi(real L, int p, vector x) {
		vector[rows(x)] fi = 1/sqrt(L) * sin(p*pi()/(2*L) * (x+L));
		return fi;
	}
}

data {
	int<lower=1> N;
	vector[N] x;
	int y[N];
	
	real L;						//boundary condition factor
	int<lower=1> P;		//number of basis functions
}

transformed data {
	matrix[N,P] PHI;
	
	for (p in 1:P) { 
	  PHI[, p] = phi(L, p, x); 
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
	
	{
	  vector[P] diagSPD;
	  vector[P] SPD_beta;
	
	  for(p in 1:P){ 
		  diagSPD[p] =  sqrt(spd(alpha, rho, sqrt(lambda(L, p)))); 
	  }
	
	  SPD_beta = diagSPD .* beta;
	
	  f = m + PHI * SPD_beta;
	}
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
  vector[N] yPred;
  vector[N] mu = exp(f);
  
  for(i in 1:N) {
    yPred[i] = f[i] > 20 ? neg_binomial_2_log_rng(20, theta)
                         : neg_binomial_2_log_rng(f[i], theta);
  }
}
