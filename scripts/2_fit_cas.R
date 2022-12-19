library(tidyverse)
library(rstan)


options(mc.cores = 4)
rstan_options(auto_write = TRUE)


m_cas_a0 <- rstan::stan_model(here::here("stan", "cas_a0.stan"))
m_cas_a <- rstan::stan_model(here::here("stan", "cas_a.stan"))
m_cas_b <- rstan::stan_model(here::here("stan", "cas_b.stan"))
m_cas_c <- rstan::stan_model(here::here("stan", "cas_c.stan"))


dir.create("out/sub_cas", showWarnings = F)


load(here::here("data", "dat_cas.rdata"))


## Exogenous variables
exo <- list(
  scale_dur = 1,
  r_death_a = 0,
  r_death_s = 0.12,
  ppv_pub = 0.85,
  ppv_eng = 0.85,
  sens_acf = 0.5,
  spec_acf = 0.995,
  dur_tx = 0.5
)


threshold <- 1.05
n_iter <- 5000
n_warmup <- 4000

locations <- unique(dat_noti$State)
locations <- locations[locations %in% c("India", "Andhra_Pradesh")]


fit <- function(m, dat, outfile, threshold=threshold) {
  rhat <- 10
  while(rhat > threshold) {
    cat(loc); cat("\t"); cat(outfile); cat("\n")
    post <- sampling(m_cas_a, data=dat, iter=n_iter, warmup=n_warmup)
    rhat <- max(summary(post)$summary[, "Rhat"], na.rm=T)
  }
  save(post, dat, file = here::here("out", "sub_cas", outfile))
  return(rhat)
}


loc <- "India"

for (loc in locations) {
  sel_noti <- dat_noti %>% filter(State == loc) %>% 
    filter(!is.na(N_Det_Pub)) %>% 
    filter(N_Det_ACF > 0)
  sel_tbps <- dat_tbps %>% filter(State == loc)
  
  
  # No drug sale data
  dat <- list(
    YearSurveyed = 2020,
    N = as.integer(sel_tbps$N),
    Asym = as.integer(sel_tbps$N_Asym),
    Sym = as.integer(sel_tbps$N_NotCS),
    CS = as.integer(sel_tbps$N_NotDet),
    Tx = as.integer(sel_tbps$N_OnATT_Pri) + as.integer(sel_tbps$N_OnATT_Pub),
    n_t = nrow(sel_noti),
    Pop = sel_noti$Pop,
    NotiPub = sel_noti$N_Det_Pub,
    NotiEng = sel_noti$N_Det_Eng,
    NotiACF = sel_noti$N_Det_ACF,
    Years = sel_noti$Year
  )
  
  dat <- c(dat, exo)
  
  fit(m_cas_a0, dat, "post_cas_a0_" + glue::as_glue(loc) + ".rdata", threshold)
  fit(m_cas_a, dat, "post_cas_a_" + glue::as_glue(loc) + ".rdata", threshold)
  fit(m_cas_b, dat, "post_cas_b_" + glue::as_glue(loc) + ".rdata", threshold)
  fit(m_cas_c, dat, "post_cas_c_" + glue::as_glue(loc) + ".rdata", threshold)
}




