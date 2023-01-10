library(tidyverse)



load(here::here("data", "dat_cas.rdata"))

locations <- unique(dat_noti$State)


## Setup
set.seed(1166)

n_samples <- 300

dir.create(here::here("docs", "pars"), showWarnings = F)

## Collect pars
for (loc in locations[c(1, 2, 5, 8)]) {
  loc <- glue::as_glue(loc)
  
  
  if (loc == "Arunachal_Pradesh") {
    txi <- jsonlite::read_json(here::here("out", "sub_txi", "pars_India.json"))
    txi <- txi[sample(1:length(txi), n_samples, rep =T)]
  } else {
    txi <- jsonlite::read_json(here::here("out", "sub_txi", "pars_" + loc + ".json"))
    txi <- txi[sample(1:length(txi), n_samples, rep =T)]
  }

  
  txo <- jsonlite::read_json(here::here("out", "sub_tx", "pars_" + loc + ".json"))
  txo <- txo[sample(1:length(txo), n_samples, rep =T)]
  

  cs <- local({
    load(here::here("out", "sub_csd", "post_csd_d2_" + loc + ".rdata"))
    
    rstan::extract(post, 
                   pars =  c(
                     "r_death_a", "r_acf", "r_sc",
                     "r_onset", "r_csi", "r_recsi", "r_acf", 
                      "ppv_pri", "adr",
                      "p_dx1_pub", "p_dx1_eng", "p_dx1_pri", 
                      "p_dx0_pub", "p_dx0_eng", "p_dx0_pri",
                      "txi_pri", "p_pri_on_pub", "ppv_pri",
                      "prv_a", "prv_s", "prv_c", "inc0",
                      "p_pub", "p_eng", "p_pri", "wt[1]"
                     )) %>% 
      as_tibble() %>% 
      rename(wt = "wt[1]", inc = inc0) %>% 
      mutate(across(c(starts_with("prv_"), inc), function(x) x * wt)) %>% 
      select(-wt) %>% 
      bind_cols(as_tibble(dat[c("r_death_s", "ppv_pub", "ppv_eng", "sens_acf", "spec_acf", "dur_pub")])) %>% 
      mutate(
        Key = 1:n()
      ) %>% 
      filter(Key %in% sample(1:n(), n_samples, rep =T)) %>% select(-Key)
  })
  
  
  pars <- lapply(1:n_samples, function(i) c(as.list(cs[i, ]), txo[[i]], txi[[i]]))
  
  jsonlite::write_json(pars, here::here("docs", "pars", "pars_" + loc + ".json"), auto_unbox=T, digits = 10)
}

