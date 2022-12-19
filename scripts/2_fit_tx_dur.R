library(tidyverse)
library(rstan)


options(mc.cores = 4)
rstan_options(auto_write = TRUE)


dir.create("out/sub_txdur", showWarnings = F)


## Models
# S1: drug sale
# S2: On ATT, TBPS
# S3: drug sale + on ATT


m_dur <- rstan::stan_model(here::here("stan", "txdur.stan"))



## Data loading
p_txi <- read_csv(here::here("docs", "tabs", "post_txi.csv"))
p_txi


load(here::here("data", "dat_cas.rdata"))

dat_noti

dat_tbps


## Exogenous variables


threshold <- 1.05
n_iter <- 5000
n_warmup <- 4000

locations <- unique(dat_noti$State)
locations


loc <- "India"

sel_noti <- dat_noti %>% filter(State == loc) %>% filter(!is.na(N_Det_Pub))
sel_txi <- p_txi %>% filter(State == loc)

dat <- list(
  n_t = nrow(sel_noti),
  Pop = sel_noti$Pop,
  NotiEng = sel_noti$N_Det_Eng,
  Years = sel_noti$Year,
  txi_eng = sel_txi$p_txi_eng_m,
  DrugPri = sel_tbps$DrugTime_M,
  DrugPri_Std = (sel_tbps$DrugTime_U - sel_tbps$DrugTime_L) / 2 / 1.96
)


post <- sampling(m_dur, data=dat, iter=n_iter, warmup=n_warmup)



while (max(rhat) > threshold) {
  for (loc in locations) {

  
    # No drug sale data
    dat <- list(
      YearSurveyed = 2020,
      N = as.integer(sel_tbps$N),
      Asym = as.integer(sel_tbps$N_Asym),
      Sym = as.integer(sel_tbps$N_NotCS),
      CS = as.integer(sel_tbps$N_NotDet),
      TxPri = as.integer(sel_tbps$N_OnATT_Pri),
      TxPub = as.integer(sel_tbps$N_OnATT_Pub),
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
    if (rhat0 > threshold) {
      cat(loc); cat("\tS1\t")
      post <- sampling(m_cs_s1, data=dat, iter=n_iter, warmup=n_warmup)
      rhat1 <- max(summary(post)$summary[, "Rhat"])
      
      if (rhat1 < rhat0) {
        rhat[loc, 1] <- rhat1
        save(post, file = here::here("out", "sub_cs", "post_cs_s1_" + glue::as_glue(loc) + ".rdata"))
      }
      cat(min(rhat1, rhat0)); cat("\n")
    }
    
    rhat0 <- rhat[loc, 2]
    if (rhat0 > threshold) {
      cat(loc); cat("\tS2\t")
      post <- sampling(m_cs_s2, data=dat, iter=n_iter, warmup=n_warmup)
      rhat1 <- max(summary(post)$summary[, "Rhat"])
      
      if (rhat1 < rhat0) {
        rhat[loc, 2] <- rhat1
        save(post, file = here::here("out", "sub_cs", "post_cs_s2_" + glue::as_glue(loc) + ".rdata"))
      }
      cat(min(rhat1, rhat0)); cat("\n")
    }
    
    rhat0 <- rhat[loc, 3]
    if (rhat0 > threshold) {
      cat(loc); cat("\tS3\t")
      post <- sampling(m_cs_s3, data=dat, iter=n_iter, warmup=n_warmup)
      rhat1 <- max(summary(post)$summary[, "Rhat"])
      
      if (rhat1 < rhat0) {
        rhat[loc, 3] <- rhat1
        save(post, file = here::here("out", "sub_cs", "post_cs_s3_" + glue::as_glue(loc) + ".rdata"))
      }
      cat(min(rhat1, rhat0)); cat("\n")
    }
  }
}



write_csv(as_tibble(rhat), file = here::here("docs", "tabs", "rhat_cs.csv"))

