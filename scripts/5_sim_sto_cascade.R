library(tidyverse)


source(here::here("R", "sim_des.R"))


each <- 5


for(cnr_year in 2019:2021) {
  folder <- glue::as_glue("cas_") + cnr_year
  
  for(scenario in c("shared_pr_asym", "shared_r_onset")) {
    file_shift <- glue::as_glue("shifting_") + scenario + ".rdata"
    load(file = here::here("out", folder, file_shift))
    locs <- names(pars_shifting)
    
    
    sims_pp <- bind_rows(lapply(locs, function(loc) {
      print(loc)
      pars <- pars_shifting[[loc]]
      ks <- pars$Rates %>% pull(Key) %>% unique()

      ks <- sample(ks, round(length(ks) / each))
      
      bind_rows(lapply(ks, function(k) {
        d <- pars$Rates %>% filter(Key == k) %>% as.list()
        
        inp <- c(list(
          Key = k,
          prev0 = d$Prev,
          r_sc = d$r_sc,
          r_die_asym = d$r_die_asym,
          r_die_sym = d$r_die_sym,
          r_onset = d$r_onset,
          r_aware = d$r_aware,
          r_csi = d$r_cs,
          r_recsi = d$r_recsi,
          p_txi = c(d$TxI_pub, d$TxI_eng, d$TxI_pri),
          r_succ = c(d$r_succ_pub, d$r_succ_eng, d$r_succ_pri),
          r_ltfu = c(d$r_ltfu_pub, d$r_ltfu_eng, d$r_ltfu_pri),
          r_die_tx = c(d$r_dead_pub, d$r_dead_eng, d$r_dead_pri)
        ), pars$Shifting)
        
        bind_rows(lapply(1:each, function(i) {
          sample_des(i, inp = inp)
        })) %>% mutate(Key = i)
      })) %>% mutate(Location = loc)
    })) %>% 
      mutate(
        SectorStart = factor(SectorStart, c("pub", "eng", "pri")),
        SectorEnd = factor(SectorEnd, c("pub", "eng", "pri"))
      )
    
    
    save(sims_pp, file = here::here("out", folder, "sims_pp_" + glue::as_glue(scenario) + ".rdata"))
  }
}

