library(tidyverse)
library(jsonlite)



p_dx1_pub = 0.7
rat_dx01 = 0.1


load(here::here("data", "dat_cas.rdata"))
dir.create(here::here("out", "pars"), showWarnings = F)


for (loc in dat_tbps$State) {
  loc <- glue::as_glue(loc)

  
  sel_cs <- read_csv(here::here("out",  "sub_cs", "post_cs_s1_" + loc + ".csv")) %>% 
    mutate(Key = sample(1:n(), n(), rep=F)) %>% 
    arrange(Key) %>% 
    select(Key, adr, prv0, 
           r_death_a, r_death_s, r_sc,
           r_sym, r_csi = r_aware,
           r_det, r_det_all, 
           det_pub, det_eng, det_pri, ppv_pub, ppv_eng, ppv_pri, p_pri_on_pub,
           r_acf, sens_acf, spec_acf,
           dur_pub, dur_pri,
           p_txi_pri=txi_pri)
  
  if (loc == "Arunachal_Pradesh") {
    sel_txi <- bind_rows(read_json(here::here("out", "sub_txi", "pars_India.json"))) %>% 
      mutate(Key = sample(1:n(), n(), rep=F)) %>% 
      arrange(Key)
  } else {
    sel_txi <- bind_rows(read_json(here::here("out", "sub_txi", "pars_" + loc + ".json"))) %>% 
      mutate(Key = sample(1:n(), n(), rep=F)) %>% 
      arrange(Key)
  }

  
  sel_txo <- bind_rows(read_json(here::here("out", "sub_tx", "pars_" + loc + ".json"))) %>% 
    mutate(Key = sample(1:n(), n(), rep=F)) %>% 
    arrange(Key)
  
  
  pars <- sel_cs %>% 
    left_join(sel_txi) %>% 
    left_join(sel_txo) %>%
    mutate(
      r_det_all = r_det * (det_pub * p_txi_pub + det_eng * p_txi_eng + det_pri * p_txi_pri),
      out_a = r_acf * sens_acf + r_death_a + r_sc + r_sym,
      out_s = r_acf * sens_acf + r_death_s + r_sc + r_csi,
      out_c = r_acf * sens_acf + r_death_s + r_sc + r_det_all,
      prv_a = 1,
      prv_s = prv_a * r_sym / (out_s - adr),
      prv_c = prv_s * r_csi / (out_c - adr),
      ps = prv_a + prv_s + prv_c,
      prv_a = prv_a / ps * prv0,
      prv_s = prv_s / ps * prv0,
      prv_c = prv_c / ps * prv0
    ) %>% 
    select(-out_a, -out_s, -out_c, -ps) %>%
    mutate(
      p_dx1_pub = p_dx1_pub, 
      p_dx1_pri = runif(n(), 0, p_dx1_pub),
      p_dx1_eng = (p_dx1_pub + p_dx1_pri) / 2,
      p_dx0_pub = p_dx1_pub * rat_dx01,
      p_dx0_eng = p_dx1_eng * rat_dx01,
      p_dx0_pri = p_dx1_pri * rat_dx01,
      # txi2 = prv_c * r_det * (det_pub * p_txi_pub + det_eng * p_txi_eng + det_pri * p_txi_pri),
      p_ent_pub = det_pub / p_dx1_pub,
      p_ent_eng = det_eng / p_dx1_eng,
      p_ent_pri = det_pri / p_dx1_pri,
      sk = p_ent_pub + p_ent_eng + p_ent_pri,
      p_ent_pub = p_ent_pub / sk,
      p_ent_eng = p_ent_eng / sk,
      p_ent_pri = p_ent_pri / sk,
      p_txi0 = (p_ent_pub * p_dx0_pub * p_txi_pub + 
                  p_ent_eng * p_dx0_eng * p_txi_eng + 
                  p_ent_pri * p_dx0_pri * p_txi_pri),
      p_txi1 = (p_ent_pub * p_dx1_pub * p_txi_pub + 
                  p_ent_eng * p_dx1_eng * p_txi_eng + 
                  p_ent_pri * p_dx1_pri * p_txi_pri),
      r_recsi = (prv_c * r_det_all - r_csi * p_txi0 * prv_s) / (p_txi1 * prv_c)
      # txi0 = r_csi * p_txi0 * prv_s,
      # txi1 = r_recsi * p_txi1 * prv_c,
      # txi_adj = txi0 + txi1,
      # prop0 = txi0 / txi_adj,
      # prop1 = txi1 / txi_adj
    ) %>% 
    select(-r_det_all, -sk)
  
  
  write_csv(pars, here::here("out", "pars", "pars_" + loc + ".csv"))
}


