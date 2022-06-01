library(tidyverse)


source(here::here("R", "ext_rate.R"))


n_mc <- 1000

ppv <- 0.65

settings <- list(
  TxI_pri = 0.8,
  r_succ_pub = 2,
  r_succ_eng = 2,
  r_succ_pri = 2,
  rr_ltfu_pri = 1.5,
  ppv_pri = 0.5
)

exo <- tibble(
  Key = 1:n_mc,
  r_die_untx = runif(n_mc, 0.14, 0.18),
  r_sc = runif(n_mc, 0.1, 0.3),
  rr_die_asym = runif(n_mc, 0, 1)
)


for (cnr_year in 2019:2021) {
  folder <- paste0("cas_", cnr_year)
  
  dir.create(here::here("out", folder), showWarnings = F)
  
  load(here::here("data", "cascade", paste0("d_cascade_", cnr_year, ".rdata")))
  
  rates <- d_cascade %>% rename(Location = State) %>% ext_rate(exo, settings)
  
  save(rates, settings, file = here::here("out", folder, "rates.rdata"))
}


settings <- list(
  TxI_pri = 0.8,
  r_succ_pub = 2,
  r_succ_eng = 2,
  r_succ_pri = 2,
  rr_ltfu_pri = 1.5,
  ppv_pub = 0.85,
  ppv_pri = 0.2
)


d_cascade %>% rename(Location = State) %>% ext_rate(exo, settings) %>% 
  group_by(Location) %>% 
  select(starts_with("DetR_"), starts_with("TxR_"), Pr_Pub_CSI) %>%
  summarise(across(everything(), mean)) %>% 
  mutate(
    DetR = DetR_pub + DetR_eng + DetR_pri,
    across(c(DetR_pub, DetR_eng, DetR_pri), function(x) x / DetR),
    TxR = TxR_pub + TxR_eng + TxR_pri,
    across(c(TxR_pub, TxR_eng, TxR_pri), function(x) x / TxR)
  ) %>% 
  head(20)

