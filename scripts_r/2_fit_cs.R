library(tidyverse)
library(rstan)


options(mc.cores = 4)
rstan_options(auto_write = TRUE)


dir.create("out/sub_cs", showWarnings = F)


m_cs <- rstan::stan_model(here::here("stan", "cs_no_drugsale.stan"))
m_csd <- rstan::stan_model(here::here("stan", "cs_with_drugsale.stan"))


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

locations <- unique(dat_noti$State)


rhat <- matrix(0, length(locations), 2)
rownames(rhat) <- locations
colnames(rhat) <- c("NoDS", "WithDS")


for (loc in locations[c(1)]) {
  print(loc)
  sel_noti <- dat_noti %>% filter(State == loc) %>% filter(!is.na(N_Det_Pub))
  sel_tbps <- dat_tbps %>% filter(State == loc)
  sel_txi <- p_txi %>% filter(State == loc)
  
  if (nrow(sel_txi) == 0) {
    sel_txi <- p_txi %>% filter(State == "India")
  }
  
  # No drug sale data
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
  
  
  po_cs_nods <- sampling(m_cs, data = dat, iter=5000, warmup=4000)
  
  rhat[loc, 1] <- max(summary(po_cs_nods)$summary[, "Rhat"])
  
  hist(rstan::extract(po_cs_nods,  "txi_pri")[[1]])
  
  save(po_cs_nods, file = here::here("out", "sub_cs", "post_cs_nods_" + glue::as_glue(loc) + ".rdata"))
  
  res_nods <- as_tibble(rstan::extract(po_cs_nods)) %>% 
    bind_cols(bind_rows(exo)) %>% 
    select(-starts_with("nr_"), - prv, - lp__)
  
  write.csv(res_nods, file=here::here("docs", "tabs", "post_cs_nods_" + glue::as_glue(loc) + ".csv"), row.names=F)
  
  
  
  
  # No drug sale data
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
    txi_eng = sel_txi$p_txi_eng_m,
    DrugPri = sel_tbps$DrugTime_M,
    DrugPri_Std = (sel_tbps$DrugTime_U - sel_tbps$DrugTime_L) / 2 / 1.96
  )
  
  dat <- c(dat, exo)
  
  
  po_cs_onds <- sampling(m_csd, data = dat, iter=5000, warmup=4000)
  
  rhat[loc, 2] <- max(summary(po_cs_onds)$summary[, "Rhat"])
  
  hist(rstan::extract(po_cs_onds,  "txi_pri")[[1]])
  
  save(po_cs_onds, file = here::here("out", "sub_cs", "post_cs_onds_" + glue::as_glue(loc) + ".rdata"))
  
  res_onds <- as_tibble(rstan::extract(po_cs_onds)) %>% 
    bind_cols(bind_rows(exo)) %>% 
    select(-starts_with("nr_"), - prv, - lp__)
  
  write.csv(res_onds, file=here::here("docs", "tabs", "post_cs_onds_" + glue::as_glue(loc) + ".csv"), row.names=F)
  
  print(rhat[loc, ])
}


write_csv(as_tibble(rhat), file = here::here("docs", "tabs", "rhat_cs.rdata"))

