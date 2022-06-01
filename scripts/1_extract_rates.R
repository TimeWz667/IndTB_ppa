library(tidyverse)


#source(here::here("R", "calc_rate.R"))


n_mc <- 1000


settings <- list(
  TxI_pri = 0.8,
  r_succ_pub = 2,
  r_succ_eng = 2,
  r_succ_pri = 1 / 0.75
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
