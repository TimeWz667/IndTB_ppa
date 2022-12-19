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


model <- odin::odin("odin/m_cas_abc.R")
cm <- model$new()




loc <- "Andhra_Pradesh"


load(here::here("out", "sub_cas", "post_cas_b_" + glue::as_glue(loc) + ".rdata"))

exo <- dat[c("r_death_a", "r_death_s", "ppv_pub", "ppv_eng", "sens_acf", "spec_acf", "dur_tx")]
pp <- rstan::extract(post, pars=names.pars) %>% 
  as_tibble() %>% 
  rename(wt = "wt[1]") %>% 
  mutate(across(c(starts_with("prv_")), function(x) x * wt)) %>% 
  select(-wt) %>% 
  bind_cols(exo)




sims <- bind_rows(lapply(1:nrow(pp), function(i) {
  p <- as.list(pp[i, ])
  
  cm$set_user(user=p)
  cm$set_user(intv = 0)
  
  ys <- cm$run(2019:2025)
  

  c(p, ys[nrow(ys), ])
  
}))



colnames(sims)


sims %>%
  summarise(
    across(c(p_dx1_pub, p_dx0_pub, p_pub, p_eng, p_pri, dur_a, dur_s, dur_c, del_pat, del_sys, del_tot), list(
      M = mean,
      L = function(x) quantile(x, 0.025),
      U = function(x) quantile(x, 0.975)
    ))
  ) %>% 
  pivot_longer(everything()) %>% 
  extract(name, c("Index", "name"), "(\\S+)_(M|L|U)") %>% 
  pivot_wider() %>% 
  mutate(M = M)



