library(tidyverse)
library(rstan)


options(mc.cores = 4)
rstan_options(auto_write = TRUE)


m_cs <- rstan::stan_model(here::here("stan", "cs_no_drugsale.stan"))


## Data loading
p_txi <- read_csv(here::here("docs", "tabs", "post_txi.csv"))

p_txi


load(here::here("data", "dat_cas.rdata"))

dat_noti

dat_tbps

## Exogenous variables
exo <- list(
  scale_dur = 1,
  r_death_a = 0,
  r_death_s = 0.12,
  ppv_pub = 0.85,
  ppv_eng = 0.85,
  sens_acf = 0.7,
  spec_acf = 0.99,
  dur_pub = 0.5
)



loc <- "India"


sel_noti <- dat_noti %>% filter(State == loc)
sel_tbps <- dat_tbps %>% filter(State == loc)
sel_txi <- p_txi %>% filter(State == loc)

dat <- list(
  YearSurveyed = 2020,
  N = as.integer(sel_tbps$N),
  Asym = as.integer(sel_tbps$N_Asym),
  Sym = as.integer(sel_tbps$N_NotCS),
  CS = as.integer(sel_tbps$N_NotDet),
  TxPri = as.integer(sel_tbps$N_OnATT_Pri),
  n_t = nrow(sel_noti),
  Pop = sel_noti$Pop,
  NotiPub = sel_noti$N_Det_Pub,
  NotiEng = sel_noti$N_Det_Eng,
  NotiACF = sel_noti$N_Det_ACF,
  Years = sel_noti$Year,
  txi_pub = sel_txi$p_txi_pub_m,
  txi_eng = sel_txi$p_txi_eng_m
)

dat <- c(dat, exo)


po <- sampling(m_cs, data = dat, iter=5000, warmup=4000)

hist(rstan::extract(po,  "txi_pri")[[1]])





