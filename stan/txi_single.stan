data {
  int<lower=0> N;
  int<lower=0> D_Pub[N];
  int<lower=0> D_Eng[N];
  int<lower=0> T_Pub[N];
  int<lower=0> T_Eng[N];
}
parameters {
  real<lower=0, upper=1> p_pub;
  real<lower=0, upper=1> p_eng;
}
model {
  for (i in 1:N) {
    target += binomial_lpmf(T_Pub[i] | D_Pub[i], p_pub);
    target += binomial_lpmf(T_Eng[i] | D_Eng[i], p_eng);
  }
}