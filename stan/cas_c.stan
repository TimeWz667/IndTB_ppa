data {
  // Data from a prevalence survey
  int<lower=0> N;
  int<lower=0> Asym;
  int<lower=0> Sym;
  int<lower=0> CS;
  int<lower=0> Tx;
  real YearSurveyed; // timing of the survey

  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  int<lower=0> Pop[n_t]; // population size
  int<lower=0> NotiPub[n_t]; // notification counts
  int<lower=0> NotiEng[n_t]; // notification counts
  int<lower=0> NotiACF[n_t]; // notification counts
  real<lower=0> Years[n_t]; // years of the notification data

  // Prior knowledge
  real<lower=0> scale_dur;

  // Exogenous variables
  real<lower=0> r_death_a;
  real<lower=0> r_death_s;
  
  real<lower=0, upper=1> ppv_pub;
  real<lower=0, upper=1> ppv_eng;
  
  real<lower=0, upper=1> sens_acf;
  real<lower=0, upper=1> spec_acf;
  
  real<lower=0> dur_tx;
}
parameters {
  real<lower=0, upper=0.05> inc0;
  real<lower=0, upper=0.1> adr;
  
  real<lower=0, upper=1> p_pub;
  real<lower=0, upper=1> p_ppm;
  
  real<lower=0, upper=1> p_dx1_pub;
  real<lower=0, upper=1> rat_dx1_pri;
  real<lower=0, upper=1> rat_dx01;
  
  real<lower=0> r_onset;
  real<lower=0> r_csi;
  real<lower=0> r_recsi;
  real<lower=0> r_acf;
  real<lower=0.1, upper=0.3> r_sc;
  
  real<lower=0.5, upper=1> txi_pri;
  real<lower=0.1, upper=ppv_eng> ppv_pri;
}
transformed parameters {
  real<lower=0> ra;
  real<lower=0> rs;
  real<lower=0> rc;
  
  real<lower=0, upper=1> p_dx1_eng;
  real<lower=0, upper=1> p_dx1_pri;
  real<lower=0, upper=1> p_dx0_pub;
  real<lower=0, upper=1> p_dx0_eng;
  real<lower=0, upper=1> p_dx0_pri;
  
  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_s;
  real<lower=0, upper=1> prv_c;
  real<lower=0, upper=1> prv0;  
  
  real<lower=0, upper=1> prv_t;
  
  real<lower=0, upper=1> p_eng;
  real<lower=0, upper=1> p_pri;

  vector<lower=0>[n_t] wt;
  vector<lower=0>[n_t] prv;
  vector<lower=0>[n_t] nr_acf;
  vector<lower=0>[n_t] nr_pub;
  vector<lower=0>[n_t] nr_eng;
  
  real p_dx0;
  real p_dx1;
  
  ra = r_sc + r_death_a + r_acf * sens_acf;
  rs = r_sc + r_death_s + r_acf * sens_acf;
  rc = r_sc + r_death_s + r_acf * sens_acf;
  
  p_dx1_pri = p_dx1_pub * rat_dx1_pri;
  p_dx0_pub = p_dx1_pub * rat_dx01;
  p_dx0_pri = p_dx1_pri * rat_dx01;

  p_dx0_eng = (p_dx0_pub + p_dx0_pri) / 2;  
  p_dx1_eng = (p_dx1_pub + p_dx1_pri) / 2;

  p_eng = (1 - p_pub) * p_ppm;
  p_pri = (1 - p_pub) * (1 - p_ppm);
  p_dx0 = p_pub * p_dx0_pub + p_eng * p_dx0_eng + p_pri * p_dx0_pri;
  p_dx1 = p_pub * p_dx1_pub + p_eng * p_dx1_eng + p_pri * p_dx1_pri;
  
  
  prv_a = inc0 / (r_onset + ra - adr);
  prv_s = r_onset * prv_a / (r_csi + rs - adr);
  prv_c = r_csi * (1 - p_dx0) * prv_s / (r_recsi * p_dx1 + rc - adr);
  prv0 = prv_a + prv_s + prv_c;

  prv_t = (prv_s * r_csi * p_pub * p_dx0_pub + prv_c * r_recsi * p_pub * p_dx1_pub) / ppv_pub;
  prv_t += (prv_s * r_csi * p_eng * p_dx0_eng + prv_c * r_recsi * p_eng * p_dx1_eng) / ppv_eng;
  prv_t += (prv_s * r_csi * p_pri * p_dx0_pri + prv_c * r_recsi * p_pri * p_dx1_pri) / ppv_pri;
  prv_t /= (1 / dur_tx - adr);
  
  
  for (i in 1:n_t) {
    wt[i] = exp(- adr * (Years[i] - YearSurveyed));
    prv[i] = wt[i] * prv0;
    
    nr_acf[i] = (prv[i] * sens_acf + (1 - prv[i]) * (1 - spec_acf)) * r_acf;
    nr_pub[i] = wt[i] * (prv_s * r_csi * p_pub * p_dx0_pub + prv_c * r_recsi * p_pub * p_dx1_pub) / ppv_pub;
    nr_eng[i] = wt[i] * (prv_s * r_csi * p_eng * p_dx0_eng + prv_c * r_recsi * p_eng * p_dx1_eng) / ppv_eng;
  }
}
model {
  inc0 ~ uniform(0, 0.05);
  r_onset ~ inv_gamma(scale_dur, scale_dur);
  r_csi ~ inv_gamma(scale_dur, scale_dur);
  r_recsi ~ inv_gamma(scale_dur, scale_dur);

  r_sc ~ uniform(0.1, 0.3);
  adr ~ uniform(0, 0.2);


  target += binomial_lpmf(Asym | N, prv_a);
  target += binomial_lpmf(Sym | N, prv_s);
  target += binomial_lpmf(CS | N, prv_c);
  target += binomial_lpmf(Tx | N, prv_t);

  for (i in 1:n_t) {
    target += poisson_lpmf(NotiACF[i] | nr_acf[i] * Pop[i]);
    target += poisson_lpmf(NotiPub[i] | nr_pub[i] * Pop[i]);
    target += poisson_lpmf(NotiEng[i] | nr_eng[i] * Pop[i]);
  }
}
generated quantities {
  vector<lower=0>[n_t] inc;
  
  for (i in 1:n_t) {
    inc[i] = inc0 * wt[i];
  }
}
