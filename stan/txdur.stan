data {
  real<lower=0> DrugPri;
  real<lower=0> DrugPri_Std;

  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  int<lower=0> Pop[n_t]; // population size
  int<lower=0> NotiEng[n_t]; // notification counts
  real<lower=0> Years[n_t]; // years of the notification data

  // Prior knowledge
  real<lower=0, upper=1> txi_eng;
}
parameters {
  real<lower=0, upper=1> det_pri;
  
  real<lower=0, upper=0.1> adr;
  real<lower=0, upper=1> p_eng;

  real<lower=0.5, upper=1> txi_pri;
  real<lower=0.04166667, upper=2> dur_pri;
  real<lower=0, upper=1> p_pri_on_pub;
}
transformed parameters {
  real<lower=0, upper=1> drug_pri;

  vector<lower=0>[n_t] nr_e;

  for (i in 1:n_t) {
    nr_e[i] = det_pri * p_eng * exp(-adr * (Years[i] - 2020));
  }
  
  drug_pri = det_pri * ((1 - p_pri_on_pub) * p_eng * txi_eng + (1 - p_eng) * txi_pri) * dur_pri;
}
model {
  p_pri_on_pub ~ beta(1.5, 3.5);
  dur_pri ~ uniform(0, 2);
  
  adr ~ uniform(-0.1, 0.1);

  target += normal_lpdf(DrugPri | drug_pri, DrugPri_Std);
  
  for (i in 1:n_t) {
    target += poisson_lpmf(NotiEng[i] | nr_e[i] * Pop[i]);
  }
}

