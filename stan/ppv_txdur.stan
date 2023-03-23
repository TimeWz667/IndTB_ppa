data {
  int<lower=0> N_Det_Pub;
  int<lower=0> N_Det_ACF;
  int<lower=0> N_Det_Eng;
  int<lower=0> N_Txi_Pub;
  int<lower=0> N_Txi_Eng;
  int<lower=0> Pop;
  real<lower=0> p_txi_pub;
  real<lower=0> p_txi_eng;
  

  real<lower=0> Drug;
  real<lower=0> Drug_Std;
  
  int<lower=0> Tx;
  int<lower=0> Tx_Pub;
  
  real<lower=0.5, upper=1> ppv_pub;

}
parameters {
  real<lower=0.05, upper=1> ppm;
  
  real<lower=0.5, upper=p_txi_eng> p_txi_pri;
  real<lower=0.2, upper=1> ppv_eng;
  real<lower=0.2, upper=ppv_eng> ppv_pri;
  
  
  real<lower=0.04166667, upper=2> dur_pri;
  real<lower=0, upper=1> p_pri_on_pub;
  
}
transformed parameters {
  real<lower=0, upper=1> det_pub;
  real<lower=0, upper=1> det_eng;
  real<lower=0, upper=1> det_pri;
  
  real<lower=0, upper=1> txi_pub;
  real<lower=0, upper=1> txi_eng;
  real<lower=0, upper=1> txi_pri;
  
  real tx_pri;
  real drug_time;
  real<lower=0, upper=1> p_pub;
  
  det_pub = (N_Det_Pub + N_Det_ACF) * 1.0 / Pop;
  det_eng =  N_Det_Eng * 1.0 / Pop;
  det_pri = det_eng * (1 - ppm) / ppm;
  
  txi_pub = N_Txi_Pub * 1.0 / Pop * ppv_pub;
  txi_eng = N_Txi_Eng * 1.0 / Pop * ppv_eng;
  txi_pri = det_pri * p_txi_pri * ppv_pri;
  
  tx_pri = N_Txi_Eng * 1.0 / Pop * (1 - p_pri_on_pub) + txi_pri / ppv_pri;
    
  drug_time = tx_pri * dur_pri;
  p_pub = txi_pub / (txi_pub + txi_eng + txi_pri);

}
model {
  ppm ~ uniform(0, 1);

  p_pri_on_pub ~ beta(1.5, 3.5);
  
  p_txi_pri ~ uniform(0.5, p_txi_eng);
    

  target += binomial_lpmf(Tx_Pub | Tx, p_pub);
  target += normal_lpdf(Drug | drug_time, Drug_Std);
  
}

generated quantities {
  real<lower=0, upper=1> p_under;
  
  p_under = det_pri / (det_pub + det_eng + det_pri);
}

