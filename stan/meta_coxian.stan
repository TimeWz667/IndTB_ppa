data {
  int nd;
  int np;
  
  int<lower=0> N_dur[nd];
  real<lower=0> M_dur[nd];
  
  int<lower=0> N_pr[np];
  int<lower=0> X_pr[np];
  real<lower=0> X_cutoff[np];
}

parameters {
  real<lower=0, upper=1> pdx0;
  real<lower=0> rate;
}

transformed parameters {
  vector<lower=0>[nd] mu_d;
  vector<lower=0, upper=1>[np] mu_p;
  
  for (i in 1:nd) {
    mu_d[i] = (1 - pdx0) / rate;
  }
  
  for (i in 1:np) {
    mu_p[i] = pdx0 + (1 - pdx0) * (1 - exp(- rate * X_cutoff[i]));
  }
}

model {
  for (i in 1:nd) {
    target += exponential_lpdf(M_dur[i] | 1 / mu_d[i]);
  }
  
  for (i in 1:np) {
    target += binomial_lpmf(X_pr[i] | N_pr[i], mu_p[i]);
  }
}

