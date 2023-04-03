data {
  int<lower=0> N_Det_Pub;
  int<lower=0> N_Det_Eng;
  int<lower=0> N_Txi_Pub;
  int<lower=0> N_Txi_Eng;
  int<lower=0> Pop;
  
  real<lower=0> Drug;
  real<lower=0> Drug_Std;
  
  int<lower=0> Tx;
  int<lower=0> Tx_Pub;
  
  real<lower=0, upper=1> ent_pub;
  real<lower=0, upper=1> pdx_pub;
  real<lower=0.5, upper=1> ppv_pub;

}
parameters {
  real<lower=0.05, upper=1> ppm;
  
  real<lower=0.5, upper=1> p_txi_pub;
  real<lower=0.5, upper=1> p_txi_eng;
  real<lower=0.5, upper=0.8> p_txi_pri;
  
  real<lower=0, upper=1> rat_pdx_pri;
  
  real<lower=0.2, upper=ppv_pub> ppv_pri;
  
  
  real<lower=0.04166667, upper=2> dur_pri;
  real<lower=0, upper=1> p_pri_on_pub;
  
  real<lower=0, upper=1> r_cs;
}
transformed parameters {
  real<lower=0, upper=1> det_pub;
  real<lower=0, upper=1> det_eng;
  real<lower=0, upper=1> det_pri;
  
  real<lower=0, upper=1> txi_pub;
  real<lower=0, upper=1> txi_eng;
  real<lower=0, upper=1> txi_pri;
  
  real<lower=0, upper=1> pdx_eng;
  real<lower=0, upper=1> pdx_pri;
  
  real<lower=0, upper=1> ppv_eng;
  
  real tx_pri;
  real drug_time;
  real<lower=0, upper=1> p_pub;
  
  ppv_eng = ppv_pri;
  
  pdx_pri = pdx_pub * rat_pdx_pri;
  pdx_eng = pdx_pri;
  
  det_pub = r_cs * ent_pub * pdx_pub;
  det_eng = r_cs * (1 - ent_pub) * ppm * pdx_eng;
  det_pri = r_cs * (1 - ent_pub) * (1 - ppm) * pdx_pri;
  
  txi_pub = det_pub * p_txi_pub;
  txi_eng = det_eng * p_txi_eng;
  txi_pri = det_pri * p_txi_pri;
  
  tx_pri = txi_eng * (1 - p_pri_on_pub) / ppv_eng + txi_pri / ppv_pri;
    
  drug_time = tx_pri * dur_pri;
  
  p_pub = txi_pub / (txi_pub + txi_eng + txi_pri);
}
model {
  r_cs~ uniform(0, 1);
  
  ppm ~ uniform(0, 1);

  p_pri_on_pub ~ beta(1.5, 3.5);
  
  p_txi_pri ~ uniform(0.5, 0.8);
    
  target += binomial_lpmf(N_Txi_Pub | N_Det_Pub, p_txi_pub);
  target += binomial_lpmf(N_Txi_Eng | N_Det_Eng, p_txi_eng);
    
  target += binomial_lpmf(N_Det_Pub | Pop, det_pub);
  target += binomial_lpmf(N_Det_Eng | Pop, det_eng);
  
  target += binomial_lpmf(Tx_Pub | Tx, p_pub);
  
  target += normal_lpdf(Drug | drug_time, Drug_Std);
  
  
}

generated quantities {
  real<lower=0, upper=1> p_under;
  real<lower=0> tp_pri_drug;
  real<lower=0> tp_pri_drug_time;
  real<lower=0> tp_pri_txi;
  
  tp_pri_drug = (txi_eng * (1 - p_pri_on_pub) + txi_pri) * Pop;
  tp_pri_txi = (txi_eng + txi_pri) * Pop;
  tp_pri_drug_time = tp_pri_drug * dur_pri;
  
  p_under = det_pri / (det_pub + det_eng + det_pri);
}

