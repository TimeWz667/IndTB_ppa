library(tidyverse)
library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)


dir.create("out/sub_csd", showWarnings = F)


## Models
# D0: without drug scale data
# S1: with drug sale


m_cs_d0 <- rstan::stan_model("stan/cs_d0.stan")
m_cs_d1 <- rstan::stan_model("stan/cs_d1.stan")
m_cs_d2 <- rstan::stan_model("stan/cs_d2.stan")


## Data loading
p_txi <- read_csv("docs/tabs/post_txi.csv")
p_txi


load("data/dat_cas.rdata")


## Exogenous variables
exo <- list(
  scale_dur = 1,
  r_death_s = 0.2,
  ppv_pub = 0.85,
  ppv_eng = 0.85,
  sens_acf = 0.5,
  spec_acf = 0.995,
  dur_pub = 0.5
)


threshold <- 1.05
n_iter <- 5000
n_warmup <- 4000

locations <- unique(dat_noti$State)
# locations <- locations[!(locations %in% c("India", "Andhra_Pradesh"))]

rhat <- matrix(2, length(locations), 3)
rownames(rhat) <- locations
colnames(rhat) <- c("NoDS", "DS", "NoDS_1PPM")


for (loc in locations[c(1, 2, 5, 8)]) {
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
    # TxPri = as.integer(sel_tbps$N_OnATT_Pri),
    # TxPub = as.integer(sel_tbps$N_OnATT_Pub),
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
  
  rhat0 <- rhat[loc, 1]
  while (rhat0 > threshold) {
    cat(loc); cat("\tD0\t")
    post <- sampling(m_cs_d0, data=dat, iter=n_iter, warmup=n_warmup)
    rhat1 <- max(summary(post)$summary[, "Rhat"], na.rm=T)

    if (rhat1 < rhat0) {
      rhat[loc, 1] <- rhat1
      rhat0 <- rhat1
      save(dat, post, file = here::here("out", "sub_csd", "post_csd_d0_" + glue::as_glue(loc) + ".rdata"))
    }
    cat(min(rhat1, rhat0)); cat("\n")
  }
  
  rhat0 <- rhat[loc, 2]
  while (rhat0 > threshold) {
    cat(loc); cat("\tD1\t")
    post <- sampling(m_cs_d1, data=dat, iter=n_iter, warmup=n_warmup)
    rhat1 <- max(summary(post)$summary[, "Rhat"], na.rm=T)
    
    if (rhat1 < rhat0) {
      rhat[loc, 2] <- rhat1
      rhat0 <- rhat1
      save(dat, post, file = here::here("out", "sub_csd", "post_csd_d1_" + glue::as_glue(loc) + ".rdata"))
    }
    cat(min(rhat1, rhat0)); cat("\n")
  }
  
  rhat0 <- rhat[loc, 3]
  while (rhat0 > threshold) {
    cat(loc); cat("\tD1\t")
    post <- sampling(m_cs_d2, data=dat, iter=n_iter, warmup=n_warmup)
    rhat1 <- max(summary(post)$summary[, "Rhat"], na.rm=T)
    
    if (rhat1 < rhat0) {
      rhat[loc, 3] <- rhat1
      rhat0 <- rhat1
      save(dat, post, file = here::here("out", "sub_csd", "post_csd_d2_" + glue::as_glue(loc) + ".rdata"))
    }
    cat(min(rhat1, rhat0)); cat("\n")
  }
}




write_csv(as_tibble(rhat), file = here::here("docs", "tabs", "rhat_csd.csv"))



## Convert to json

for (loc in locations) {
  load(file = here::here("out", "sub_csd", "post_cs_d0_" + glue::as_glue(loc) + ".rdata"))
  
  res <- as_tibble(rstan::extract(post)) %>%
    bind_cols(bind_rows(exo)) %>%
    select(-starts_with("nr_"), - prv, - lp__)
  
  write.csv(res, file=here::here("out", "sub_csd", "post_csd_d0_" + glue::as_glue(loc) + ".csv"), row.names=F)
  
  jsonlite::write_json(apply(as.matrix(res), 1, as.list),
                       here::here("out", "sub_csd", "post_csd_d0_" + glue::as_glue(loc) + ".json"),
                       simplifyVector=T, auto_unbox=T, digits=10)
  
  
  load(file = here::here("out", "sub_csd", "post_csd_d1_" + glue::as_glue(loc) + ".rdata"))
  
  res <- as_tibble(rstan::extract(post)) %>%
    bind_cols(bind_rows(exo)) %>%
    select(-starts_with("nr_"), - prv, - lp__)
  
  write.csv(res, file=here::here("out", "sub_csd", "post_csd_d1_" + glue::as_glue(loc) + ".csv"), row.names=F)
  
  jsonlite::write_json(apply(as.matrix(res), 1, as.list),
                       here::here("out", "sub_csd", "post_csd_d1_" + glue::as_glue(loc) + ".json"),
                       simplifyVector=T, auto_unbox=T, digits=10)
  
}
