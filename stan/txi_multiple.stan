data {
  int<lower=0> N;
  int<lower=0> D_Pub[N];
  int<lower=0> D_Eng[N];
  int<lower=0> T_Pub[N];
  int<lower=0> T_Eng[N];
  
  int States[N];
  int Regions[N];
  
  int<lower=0> N_State;
  int<lower=0> N_Region;
  
  int States_Pred[N_State];
  int Regions_Pred[N_State];
}
parameters {
  real alpha;
  real alpha_state[N_State];
  real alpha_region[N_Region];
  
  real beta;
  real beta_state[N_State];
  real beta_region[N_Region];
}
model {
  for (i in 1:N) {
    target += binomial_logit_lpmf(T_Pub[i] | D_Pub[i], alpha + alpha_state[States[i]] + alpha_region[Regions[i]]);
    target += binomial_logit_lpmf(T_Eng[i] | D_Eng[i], beta + beta_state[States[i]] + beta_region[Regions[i]]);
  }
}
generated quantities {
  real<lower=0, upper=1> p_pub[N_State];
  real<lower=0, upper=1> p_eng[N_State];
  
  for (i in 1:N_State) {
    p_pub[i] = inv_logit(alpha + alpha_state[States_Pred[i]] + alpha_region[Regions_Pred[i]]);
    p_eng[i] = inv_logit(beta + beta_state[States_Pred[i]] + beta_region[Regions_Pred[i]]);
  }
}
