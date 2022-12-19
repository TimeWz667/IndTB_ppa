library(odin)
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


n_sel <- 1000
i_sel <- sort(sample(4000, n_sel))



model <- odin::odin("odin/m_cas_abc.R", target = "r")
cm <- model$new()


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


for (loc in c("India", "Andhra_Pradesh")) {
  pars <- bind_rows(
    load_pars("post_cas_a0_", loc, "Juan"),
    load_pars("post_cas_a_", loc, "Eq"),
    load_pars("post_cas_b_", loc, "AllRe"),
    load_pars("post_cas_c_", loc, "Free")
  ) %>% 
    filter(Key %in% i_sel)
  
  
  sims <- bind_rows(lapply(1:nrow(pars), function(i) {
    p <- as.list(pars[i, ])
    
    cm$set_user(user=p)
    cm$set_user(intv = 0)
    
    ys <- cm$run(2019:2025)
    
    bind_cols(p, ys)
  }))
  
  write_csv(sims, here::here("out", "sub_cas", "sims_baseline_" + glue::as_glue(loc) + ".csv"))
}



