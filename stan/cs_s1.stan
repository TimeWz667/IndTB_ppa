data {
  // Data from a prevalence survey
  int<lower=0> N;
  int<lower=0> Asym;
  int<lower=0> Sym;
  int<lower=0> CS;
  int<lower=0> TxPub;
  int<lower=0> TxPri;
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
  real<lower=0, upper=1> prv0;
  real<lower=0, upper=0.1> adr;
  real<lower=0, upper=1> p_pub;
  real<lower=0, upper=1> p_ppm;
  real<lower=0> r_sym;
  real<lower=0> r_aware;
  real<lower=0> r_acf;
  real<lower=0> r_det_all;
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
  
  real<lower=0> r_det;

  real<lower=0> a0;
  real<lower=0> c0;
  
  real<lower=0, upper=1> pr_a;
  real<lower=0, upper=1> pr_s;
  real<lower=0, upper=1> pr_c;
  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_s;
  real<lower=0, upper=1> prv_c;
  
  real<lower=0, upper=1> prv_t_pub;
  real<lower=0, upper=1> prv_t_eu;
  real<lower=0, upper=1> prv_t_ei;
  real<lower=0, upper=1> prv_t_pri;
  
  real<lower=0, upper=1> tbps_pub;
  real<lower=0, upper=1> tbps_pri;
  real<lower=0, upper=1> drug_pri;
  
  real<lower=0, upper=1> det_pub;
  real<lower=0, upper=1> det_eng;
  real<lower=0, upper=1> det_pri;

  vector<lower=0>[n_t] prv;
  vector<lower=0>[n_t] nr_acf;
  vector<lower=0>[n_t] nr_pub;
  vector<lower=0>[n_t] nr_eng;
  
  
  ra = r_sc + r_death_a + r_acf * sens_acf;
  rs = r_sc + r_death_s + r_acf * sens_acf;
  rc = r_sc + r_death_s + r_acf * sens_acf;
  
  a0 = (rs + r_aware - adr) / r_sym;
  c0 = r_aware / (rc + r_det_all - adr);
    
  pr_a = a0 / (a0 + 1 + c0);
  pr_s = 1 / (a0 + 1 + c0);
  pr_c = c0 / (a0 + 1 + c0);
  
  prv_a = prv0 * pr_a;
  prv_s = prv0 * pr_s;
  prv_c = prv0 * pr_c;
  
  det_pub = p_pub;
  det_eng = (1 - p_pub) * p_ppm;
  det_pri = (1 - p_pub) * (1 - p_ppm);
  
  r_det = r_det_all / (det_pub * txi_pub + det_eng * txi_eng + det_pri * txi_pri);
  
  for (i in 1:n_t) {
    prv[i] = prv0 * exp(- adr * (Years[i] - YearSurveyed));
    
    nr_acf[i] = (prv[i] * sens_acf + (1 - prv[i]) * (1 - spec_acf)) * r_acf;
    nr_pub[i] = prv[i] * pr_c * r_det * det_pub / ppv_pub;
    nr_eng[i] = prv[i] * pr_c * r_det * det_eng / ppv_eng;
    
  }
  
  prv_t_pub = prv0 * pr_c * r_det * det_pub * txi_pub / ppv_pub;
  prv_t_pub += (prv0 * sens_acf + (1 - prv0) * (1 - spec_acf)) * r_acf;
  prv_t_pub /= (1 / dur_pub - adr);
  prv_t_eu = prv0 * pr_c * r_det * det_eng * txi_eng / ppv_eng * p_pri_on_pub * (1 / dur_pub - adr);
  prv_t_ei = prv0 * pr_c * r_det * det_eng * txi_eng / ppv_eng * (1 - p_pri_on_pub) * (1 / dur_pri - adr);
  prv_t_pri = prv0 * pr_c * r_det * det_pri * txi_pri / ppv_pri * (1 / dur_pri - adr);
  
  tbps_pub = prv_t_pub;
  tbps_pri = (prv_t_eu + prv_t_ei + prv_t_pri);
  drug_pri = (prv_t_ei + prv_t_pri);
}
model {
  prv0 ~ uniform(0, 1);
  r_det_all ~ inv_gamma(scale_dur, scale_dur);

  r_sym ~ inv_gamma(scale_dur, scale_dur);
  r_aware ~ inv_gamma(scale_dur, scale_dur);
  r_sc ~ uniform(0.1, 0.3);

  p_pri_on_pub ~ beta(1.5, 3.5);
  
  adr ~ uniform(-0.1, 0.1);


  target += binomial_lpmf(Asym | N, prv_a);
  target += binomial_lpmf(Sym | N, prv_s);
  target += binomial_lpmf(CS | N, prv_c);
  
  // target += binomial_lpmf(TxPub | N, tbps_pub);
  // target += binomial_lpmf(TxPri | N, tbps_pri);
  target += normal_lpdf(DrugPri | (prv_t_ei + prv_t_pri), DrugPri_Std);
  
  for (i in 1:n_t) {
    target += poisson_lpmf(NotiACF[i] | nr_acf[i] * Pop[i]);
    target += poisson_lpmf(NotiPub[i] | nr_pub[i] * Pop[i]);
    target += poisson_lpmf(NotiEng[i] | nr_eng[i] * Pop[i]);
  }
}
