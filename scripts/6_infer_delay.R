library(tidyverse)
library(rstan)

names.pars <- c(
  "r_onset", "r_csi", "r_recsi", "r_acf", "r_sc", 
  "ppv_pri", "adr",
  "p_dx1_pub", "p_dx1_eng", "p_dx1_pri", 
  "p_dx0_pub", "p_dx0_eng", "p_dx0_pri", 
  "prv_a", "prv_s", "prv_c", "prv_t", 
  "p_pub", "p_eng", "p_pri", "wt[1]"
)


load_pars <- function(m, loc, label) {
  load(here::here("out", "sub_cas", m + glue::as_glue(loc) + ".rdata"))
  
  exo <- dat[c("r_death_a", "r_death_s", "ppv_pub", "ppv_eng", 
               "sens_acf", "spec_acf", "dur_tx")]
  
  pp <- rstan::extract(post, pars=names.pars) %>% 
    as_tibble() %>% 
    rename(wt = "wt[1]") %>% 
    mutate(across(c(starts_with("prv_")), function(x) x * wt)) %>% 
    select(-wt) %>% 
    bind_cols(exo) %>% 
    mutate(
      Key = 1:n(),
      Model = label
    )
  
  pp
}




for (loc in c("Andhra_Pradesh", "India")) {
  ps <- load_pars("post_cas_a0_", loc, "Eq0")
  
  tab <- ps %>% 
    mutate(
      p_dx0 = p_dx0_pub * p_pub + p_dx0_eng * p_eng + p_dx0_pri * p_pri,
      p_dx1 = p_dx1_pub * p_pub + p_dx1_eng * p_eng + p_dx1_pri * p_pri,
      
      ppm_cs = p_eng / (p_eng + p_pri),
      ppm_det = p_dx1_eng * p_eng / (p_dx1_eng * p_eng + p_dx1_pri * p_pri),
      
      mu_a = r_sc + r_death_a,
      mu_s = r_sc + r_death_s,
      mu_c = mu_s,
      
      out_a = mu_a + r_onset + r_acf * sens_acf,
      out_s = mu_s + r_csi + r_acf * sens_acf,
      out_c = mu_c + r_recsi * p_dx1 + r_acf * sens_acf,
      
      dur_a = 1 / out_a * 12,
      dur_s = 1 / out_s * 12,
      dur_c = 1 / out_c * 12,
      
      p_prog_a = r_onset / out_a,
      p_acf_a = r_acf * sens_acf / out_a,
      p_pcf_a = 0,
      p_det_a = p_acf_a + p_pcf_a,
      p_drop_a = 1 - p_prog_a - p_det_a,
      
      p_prog_s = r_csi * (1 - p_dx0) / out_s,
      p_acf_s = r_acf * sens_acf / out_s,
      p_pcf_s = r_csi * p_dx0 / out_s,
      p_det_s = p_acf_s + p_pcf_s,
      p_drop_s = 1 - p_prog_s - p_det_s,
      
      p_acf_c = r_acf * sens_acf / out_c,
      p_pcf_c = r_recsi * p_dx1 / out_c,
      p_det_c = p_acf_c + p_pcf_c,
      p_drop_c = 1 - p_det_c,
      
      cas_det_a = p_det_a,
      cas_acf_a = p_acf_a,
      cas_pcf_a = p_pcf_a,
      cas_drop_a = p_drop_a,
      
      cas_det_s = p_prog_a * p_det_s,
      cas_acf_s = p_prog_a * p_acf_s,
      cas_pcf_s = p_prog_a * p_pcf_s,
      cas_drop_s = p_prog_a * p_drop_s,
      
      cas_det_c = p_prog_a * p_prog_s * p_det_c,
      cas_acf_c = p_prog_a * p_prog_s * p_acf_c,
      cas_pcf_c = p_prog_a * p_prog_s * p_pcf_c,
      cas_drop_c = p_prog_a * p_prog_s * p_drop_c,
      
      # cas_all = cas_det_a + cas_drop_a + cas_det_s + cas_drop_s + cas_det_c + cas_drop_c,
      
      p_at0 = cas_det_s / (cas_det_s + cas_det_c),
      p_at1 = 1 - p_at0,
      
      delay_pat = dur_s,
      
      delay_sys = 0 * p_at0 + dur_c * p_at1,
      
      vis = 1 + (1 / p_dx1) * p_at1,
      vis_pub = 1 + (1 / p_dx1_pub) * p_at1
    ) %>% 
    select(#starts_with("mu_"), 
      p_dx0, p_dx0_pub, vis, vis_pub, p_at0, p_at1,
      dur_a, dur_s, dur_c,
      starts_with("p_drop_"),
      starts_with("p_det_"),
      starts_with("p_prog_"),
      starts_with("delay_"),
      starts_with("cas_"),
      starts_with("ppm_")
    ) %>% 
    pivot_longer(everything()) %>% 
    group_by(name) %>% 
    summarise_all(list(
      M = median,
      L = function(x) quantile(x, 0.25),
      U = function(x) quantile(x, 0.75) 
    ))
  
  
  tab %>% 
    mutate(
      M = ifelse(startsWith(name, "cas_") | startsWith(name, "p_") | startsWith(name, "ppm_"), scales::percent(M, 0.01),scales::number(M, accuracy = 0.1)),
      L = ifelse(startsWith(name, "cas_") | startsWith(name, "p_") | startsWith(name, "ppm_"), scales::percent(L, 0.01),scales::number(L, accuracy = 0.1)),
      U = ifelse(startsWith(name, "cas_") | startsWith(name, "p_") | startsWith(name, "ppm_"), scales::percent(U, 0.01),scales::number(U, accuracy = 0.1)),
      MLU = paste0(M, " (", L, " - ", U, ")")
    ) %>% 
    write_csv(here::here("out", "sub_cas", "tab_" + glue::as_glue(loc) + ".csv"))
}






k <- runif(5e5)
r <- rgamma(5e5, 0.1, 0.1)
nv <- rgeom(length(k), k) + 1
del <- rexp(length(k), k * r)


sel <- (nv >= 0 & nv <= 1) & (del >= 1 & del <= 2)


hist(k[sel], freq = F)
hist(r[sel], freq = F, xlim = c(0, 5), breaks = 1000)


lines(seq(0, 5, 0.1), dexp(seq(0, 5, 0.1), 0.9))



