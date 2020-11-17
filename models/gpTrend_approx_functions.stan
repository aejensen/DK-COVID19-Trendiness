/* 
  Latent Gaussian Process regression model with negative binomial outcome model
  using the Hilbert space reduced rank method of Solin and Särkkä (2020) and 
  Riutort-Mayol, Bürkner, Andersen, Solin and Vehtari (2020). 
   
  Andreas Kryger Jensen, 2020. 
*/

functions {
  real besselI(real x, real nu) {
    real half_x;
    real log_half_x;
    real initial;
    real summand;
    real biggest;
    real smallest;
    real piece;
    int m;
    int mp1;
    
    half_x = 0.5 * x;
    log_half_x = log(half_x);
    m = 0;
    mp1 = 1;
    initial = 0;
    while ((mp1 + nu) < 0) {
      initial = initial + half_x ^ (2 * m + nu) / (tgamma(mp1) * tgamma(mp1 + nu));
      m = mp1;
      mp1 = m + 1;
    }
    biggest = -lgamma(mp1) -lgamma(mp1 + nu) + (2 * m + nu) * log_half_x;
    m = mp1;
    piece = positive_infinity();
    smallest = -745.13321911;
    summand = 0.0;
    while (piece > smallest) {
      mp1 = m + 1;
      piece = -lgamma(mp1) - lgamma(mp1 + nu) + (2 * m + nu) * log_half_x - biggest;
      summand = summand + exp(piece);
      m = mp1;
    }
    return exp(biggest + log1p(summand)) + initial;
  }
  
  real besselK(real x, real nu) {
    return 0.5 * pi() * (besselI(x, -nu) - besselI(x, nu)) / sin(nu * pi());
  }
  
  real spd_SE(real alpha, real rho, real w) {
    /* Spectral density function for SE covariance function */
    return square(alpha) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
  }
  
  real spd_Matern(real alpha, real rho, real nu, real w) {
  	real k1 = (nu + 1)*log(2) + log(sqrt(pi())*square(alpha)) + nu*log(2) - 2*nu*log(rho);
  	real k2 = (-0.5-nu)*log((2*nu)/(square(rho)) + 4*square(pi() * w));
  	real logSPD = k1 + k2 + lgamma(0.5 + nu) - lgamma(nu);
  	
  	return exp(logSPD);
  }
  
  real spd_Matern52(real alpha, real rho, real w) {
  	real logSPD = 2*log(alpha) + log(400*sqrt(5)) + log(rho) - log(3) - 3*log(5 + 4*square(pi()*rho*w));
  	return exp(logSPD);
  }
  
  real spd_Matern32(real alpha, real rho, real w) {
  	real logSPD = 2*log(alpha) + log(12*sqrt(3)*rho) - 2*log(3 + 4*square(pi()*rho*w));
  	return exp(logSPD);
  }

  real spd_RQ(real alpha, real rho, real nu, real w) {
  	/* Not tested */
  	real k1 = 0.25*(7-2*nu)*log(2) + log(sqrt(pi())*square(alpha)) + 
              0.25*(1+2*nu)*log(nu*square(rho)) + 0.25*(2*nu-1)*log(square(w));
  	real k2 = lgamma(nu);
  	real k3 = log(besselK(sqrt(2)*rho*fabs(w) / sqrt(1/nu), -0.5 + nu));
  	return exp(k1 - k2 + k3);
  }
  
  real lambda(real L, int p) {
    /* Eigenvalues of Laplacian */
    real lam = ((p*pi()) / (2*L))^2;
    return lam;
  }

  vector phi(real L, int p, vector x) {
    /* Eigenvectors of Laplacian */
    vector[rows(x)] fi = 1/sqrt(L) * sin(p*pi() / (2*L) * (x + L));
    return fi;
  }
  
  vector phi_d1(real L, int p, vector x) {
    /* First-order derivatives of eigenvectors of Laplacian */
    vector[rows(x)] fi = (p*pi())/(2*pow(L, 1.5)) * cos(p*pi() / (2*L) * (x + L));
    return fi;
  }
  
  vector phi_d2(real L, int p, vector x) {
    /* Second-order derivatives of eigenvectors of Laplacian */
    vector[rows(x)] fi = -(square(p)*square(pi()))/(4*pow(L, 2.5)) * sin(p*pi() / (2*L) * (x + L));
    return fi;
  }
}
