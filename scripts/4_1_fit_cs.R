library(tidyverse)
library(rstan)


options(mc.cores = 4)
rstan_options(auto_write = TRUE)


dir.create("out/sub_cs", showWarnings = F)


## Models
# S1: drug sale
# S2: On ATT, TBPS
# S3: drug sale + on ATT

m_cs_s0 <- rstan::stan_model(here::here("stan", "cs_s0.stan"))
m_cs_s1 <- rstan::stan_model(here::here("stan", "cs_s1.stan"))
m_cs_s2 <- rstan::stan_model(here::here("stan", "cs_s2.stan"))
m_cs_s3 <- rstan::stan_model(here::here("stan", "cs_s3.stan"))


## Data loading
p_txi <- read_csv(here::here("docs", "tabs", "post_txi.csv"))
p_txi


load(here::here("data", "dat_cas.rdata"))

dat_noti

dat_tbps



dat_noti %>% 
  left_join(dat_tbps) %>% 
  filter(Year == 2020) %>% 
  mutate(
    Dur = DrugTime_M / (N_Txi_Eng / Pop * 0.7),
    PrevTxPubITR = (N_Txi_Pub + N_Txi_Eng) / Pop * 0.5,
    PrevTxPubTBPS = (N_OnATT_Pub + N_OnATT_Pri) / N
  ) %>% select(State, Dur) %>% 
  data.frame()


dat_noti %>% 
  left_join(dat_tbps) %>% 
  filter(Year == 2020) %>% 
  mutate(
    Dur = DrugTime_M / (N_Txi_Eng / Pop * 0.7),
  ) %>% select(State, Dur) %>% 
  data.frame()



dat_noti %>% 
  left_join(dat_tbps) %>% 
  filter(Year == 2020) %>% 
  mutate(
    PrevTxPubITR = N_Txi_Eng / Pop * 0.5,
    PrevTxPubTBPS = N_OnATT_Pri / N
  ) %>% 
  select(State, PrevTxPubITR, PrevTxPubTBPS, starts_with("DrugTime")) %>% 
  ggplot() +
  geom_point(aes(x = PrevTxPubITR, y = PrevTxPubTBPS)) +
  geom_text(aes(x = PrevTxPubITR, y = PrevTxPubTBPS, label = State)) + 
  geom_abline(slope = 1) + 
  expand_limits(x = 0, y = 0)



dat_noti %>% 
  left_join(dat_tbps) %>% 
  filter(Year == 2020) %>% 
  mutate(
    PrevTxPubITR = N_Txi_Pub / Pop * 0.5,
    PrevTxPubTBPS = N_OnATT_Pub / N,
    PrevTxPriDiff = N_OnATT_Pri / N + (PrevTxPubTBPS - PrevTxPubITR)
  ) %>% 
  select(State, PrevTxPubITR, PrevTxPubTBPS, PrevTxPriDiff, starts_with("DrugTime")) %>% 
  ggplot() +
  geom_point(aes(x = PrevTxPriDiff, y = DrugTime_L)) +
  geom_text(aes(x = PrevTxPriDiff, y = DrugTime_L, label = State)) + 
  geom_abline(slope = 1) + 
  expand_limits(x = 0, y = 0)


## Exogenous variables
exo <- list(
  scale_dur = 1,
  r_death_a = 0,
  r_death_s = 0.12,
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
locations

rhat <- matrix(2, length(locations), 4)
rownames(rhat) <- locations
colnames(rhat) <- c("S0", "S1", "S2", "S3")


while (max(rhat) > threshold) {
  for (loc in locations) {
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
      cat(loc); cat("\tS0\t")
      post <- sampling(m_cs_s0, data=dat, iter=n_iter, warmup=n_warmup)
      rhat1 <- max(summary(post)$summary[, "Rhat"])
      
      if (rhat1 < rhat0) {
        rhat[loc, 1] <- rhat1
        save(post, file = here::here("out", "sub_cs", "post_cs_s0_" + glue::as_glue(loc) + ".rdata"))
      }
      cat(min(rhat1, rhat0)); cat("\n")
    }
    
    rhat0 <- rhat[loc, 2]
    if (rhat0 > threshold) {
      cat(loc); cat("\tS1\t")
      post <- sampling(m_cs_s1, data=dat, iter=n_iter, warmup=n_warmup)
      rhat1 <- max(summary(post)$summary[, "Rhat"])
      
      if (rhat1 < rhat0) {
        rhat[loc, 2] <- rhat1
        save(post, file = here::here("out", "sub_cs", "post_cs_s1_" + glue::as_glue(loc) + ".rdata"))
      }
      cat(min(rhat1, rhat0)); cat("\n")
    }
    
    rhat0 <- rhat[loc, 3]
    if (rhat0 > threshold) {
      cat(loc); cat("\tS2\t")
      post <- sampling(m_cs_s2, data=dat, iter=n_iter, warmup=n_warmup)
      rhat1 <- max(summary(post)$summary[, "Rhat"])

      if (rhat1 < rhat0) {
        rhat[loc, 3] <- rhat1
        save(post, file = here::here("out", "sub_cs", "post_cs_s2_" + glue::as_glue(loc) + ".rdata"))
      }
      cat(min(rhat1, rhat0)); cat("\n")
    }

    rhat0 <- rhat[loc, 4]
    if (rhat0 > threshold) {
      cat(loc); cat("\tS3\t")
      post <- sampling(m_cs_s3, data=dat, iter=n_iter, warmup=n_warmup)
      rhat1 <- max(summary(post)$summary[, "Rhat"])

      if (rhat1 < rhat0) {
        rhat[loc, 4] <- rhat1
        save(post, file = here::here("out", "sub_cs", "post_cs_s3_" + glue::as_glue(loc) + ".rdata"))
      }
      cat(min(rhat1, rhat0)); cat("\n")
    }
  }
}



write_csv(as_tibble(rhat), file = here::here("docs", "tabs", "rhat_cs.csv"))



## Convert to json

for (loc in locations) {
  load(file = here::here("out", "sub_cs", "post_cs_s0_" + glue::as_glue(loc) + ".rdata"))
  
  res <- as_tibble(rstan::extract(post)) %>%
    bind_cols(bind_rows(exo)) %>%
    select(-starts_with("nr_"), - prv, - lp__)
  
  write.csv(res, file=here::here("out", "sub_cs", "post_cs_s0_" + glue::as_glue(loc) + ".csv"), row.names=F)
  
  jsonlite::write_json(apply(as.matrix(res), 1, as.list),
                       here::here("out", "sub_cs", "post_cs_s0_" + glue::as_glue(loc) + ".json"),
                       simplifyVector=T, auto_unbox=T, digits=10)
  
  
  load(file = here::here("out", "sub_cs", "post_cs_s1_" + glue::as_glue(loc) + ".rdata"))

  res <- as_tibble(rstan::extract(post)) %>%
    bind_cols(bind_rows(exo)) %>%
    select(-starts_with("nr_"), - prv, - lp__)

  write.csv(res, file=here::here("out", "sub_cs", "post_cs_s1_" + glue::as_glue(loc) + ".csv"), row.names=F)

  jsonlite::write_json(apply(as.matrix(res), 1, as.list),
                       here::here("out", "sub_cs", "post_cs_s1_" + glue::as_glue(loc) + ".json"),
                       simplifyVector=T, auto_unbox=T, digits=10)
  
  
  load(file = here::here("out", "sub_cs", "post_cs_s2_" + glue::as_glue(loc) + ".rdata"))

  res <- as_tibble(rstan::extract(post)) %>%
    bind_cols(bind_rows(exo)) %>%
    select(-starts_with("nr_"), - prv, - lp__)

  write.csv(res, file=here::here("out", "sub_cs", "post_cs_s2_" + glue::as_glue(loc) + ".csv"), row.names=F)

  jsonlite::write_json(apply(as.matrix(res), 1, as.list),
                       here::here("out", "sub_cs", "post_cs_s2_" + glue::as_glue(loc) + ".json"),
                       simplifyVector=T, auto_unbox=T, digits=10)

  load(file = here::here("out", "sub_cs", "post_cs_s3_" + glue::as_glue(loc) + ".rdata"))

  res <- as_tibble(rstan::extract(post)) %>%
    bind_cols(bind_rows(exo)) %>%
    select(-starts_with("nr_"), - prv, - lp__)

  write.csv(res, file=here::here("out", "sub_cs", "post_cs_s3_" + glue::as_glue(loc) + ".csv"), row.names=F)

  jsonlite::write_json(apply(as.matrix(res), 1, as.list),
                       here::here("out", "sub_cs", "post_cs_s3_" + glue::as_glue(loc) + ".json"),
                       simplifyVector=T, auto_unbox=T, digits=10)
}
