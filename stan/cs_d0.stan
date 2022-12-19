data {
  // Data from a prevalence survey
  int<lower=0> N;
  int<lower=0> Asym;
  int<lower=0> Sym;
  int<lower=0> CS;
  real YearSurveyed; // timing of the survey

  real<lower=0> DrugPri;
  real<lower=0> DrugPri_Std;

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
  
  real<lower=0, upper=1> txi_pub;
  real<lower=0, upper=1> txi_eng;
  
  real<lower=0> dur_pub;
}
parameters {
  real<lower=0, upper=0.05> inc0;
  real<lower=0, upper=0.1> adr;
  real<lower=0, upper=1> p_pub;
  real<lower=0, upper=1> p_ppm;
  
  real<lower=0> r_onset;
  real<lower=0> r_csi;
  real<lower=0> r_recsi;
  real<lower=0> r_acf;
  real<lower=0.1, upper = 0.3> r_sc;
  
  real<lower=0.5, upper=1> txi_pri;
  real<lower=0.1, upper=ppv_eng> ppv_pri;
  real<lower=0.04166667, upper=1.5> dur_pri;
  real<lower=0, upper=1> p_pri_on_pub;
}
transformed parameters {
  real<lower=0> ra;
  real<lower=0> rs;
  real<lower=0> rc;
  
  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_s;
  real<lower=0, upper=1> prv_c;
  
  real<lower=0, upper=1> prv0;
  
  real<lower=0, upper=1> prv_t_pub;
  real<lower=0, upper=1> prv_t_eu;
  real<lower=0, upper=1> prv_t_ei;
  real<lower=0, upper=1> prv_t_pri;
  
  real<lower=0, upper=1> p_eng;
  real<lower=0, upper=1> p_pri;
  
  real<lower=0, upper=1> drug_pri;
  
  real<lower=0, upper=1> det_pub;
  real<lower=0, upper=1> det_eng;
  real<lower=0, upper=1> det_pri;

  vector<lower=0>[n_t] wt;
  vector<lower=0>[n_t] prv;
  vector<lower=0>[n_t] nr_acf;
  vector<lower=0>[n_t] nr_pub;
  vector<lower=0>[n_t] nr_eng;
  
  
  ra = r_sc + r_death_a + r_acf * sens_acf;
  rs = r_sc + r_death_s + r_acf * sens_acf;
  rc = r_sc + r_death_s + r_acf * sens_acf;
  

  p_eng = (1 - p_pub) * p_ppm;
  p_pri = (1 - p_pub) * (1 - p_ppm);
  
  prv_a = inc0 / (r_onset + ra - adr);
  prv_s = r_onset * prv_a / (r_csi + rs - adr);
  
  
  det_pub = r_recsi * p_pub * txi_pub;
  det_eng = r_recsi * p_eng * txi_eng;
  det_pri = r_recsi * p_pri * txi_pri;
  
  prv_c = r_csi * prv_s / (det_pub + det_eng + det_pri + rc - adr);
  prv0 = prv_a + prv_s + prv_c;
  
  for (i in 1:n_t) {
    wt[i] = exp(- adr * (Years[i] - YearSurveyed));
    
    prv[i] = prv0 * wt[i];
    
    nr_acf[i] = (prv[i] * sens_acf + (1 - prv[i]) * (1 - spec_acf)) * r_acf;
    nr_pub[i] = prv_c * r_recsi * p_pub / ppv_pub;
    nr_eng[i] = prv_c * r_recsi * p_eng / ppv_eng;
    
  }
  
  prv_t_pub = prv_c * det_pub / ppv_pub;
  prv_t_pub += (prv0 * sens_acf + (1 - prv0) * (1 - spec_acf)) * r_acf;
  prv_t_pub /= (1 / dur_pub - adr);
  prv_t_eu = prv_c * det_eng / ppv_eng * p_pri_on_pub * (1 / dur_pub - adr);
  prv_t_ei = prv_c * det_eng / ppv_eng * (1 - p_pri_on_pub) * (1 / dur_pri - adr);
  prv_t_pri = prv_c * det_pri / ppv_pri * (1 / dur_pri - adr);
  
  drug_pri = (prv_t_ei + prv_t_pri);
}
model {
  inc0 ~ uniform(0, 0.05);
  
  r_onset ~ inv_gamma(scale_dur, scale_dur);

  r_csi ~ inv_gamma(scale_dur, scale_dur);
  r_recsi ~ inv_gamma(scale_dur, scale_dur);
  
  r_acf ~ inv_gamma(scale_dur, scale_dur);

  r_sc ~ uniform(0.1, 0.3);

  p_pri_on_pub ~ beta(1.5, 3.5);
  
  adr ~ uniform(0, 0.2);


  target += binomial_lpmf(Asym | N, prv_a);
  target += binomial_lpmf(Sym | N, prv_s);
  target += binomial_lpmf(CS | N, prv_c);
  
  // target += binomial_lpmf(TxPub | N, tbps_pub);
  // target += binomial_lpmf(TxPri | N, tbps_pri);
  // target += normal_lpdf(DrugPri | (prv_t_ei + prv_t_pri), DrugPri_Std);
  
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
