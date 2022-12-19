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


transform_by_ratio <- function(p, rat_dx=0.5) {
  p <- unlist(p)
  p1 <- p
  
  ent <- p[c("p_pub", "p_eng", "p_pri")]
  pdx0 <- p[c("p_dx0_pub", "p_dx0_eng", "p_dx0_pri")]
  pdx1 <- p[c("p_dx1_pub", "p_dx1_eng", "p_dx1_pri")]
  
  det0 <- p['r_csi'] * p["prv_s"] * sum(ent * pdx0)
  fn0 <- p['r_csi'] * p["prv_s"] - det0
  det1 <- p['r_recsi'] * p["prv_c"] * sum(ent * pdx1)
  det_all <- det0 + det1
  
  det0_new <- p['r_csi'] * p["prv_s"] * sum(ent * pdx1) * rat_dx
  
  r_recsi_new <- p['r_recsi'] * (det_all - det0_new) / det1
  
  
  p1['r_recsi'] <- r_recsi_new
  p1[c("p_dx0_pub", "p_dx0_eng", "p_dx0_pri")] <- pdx1 * rat_dx
  
  return(as.list(p1)) 
}




for (loc in c("India", "Andhra_Pradesh")) {
  load(here::here("out", "sub_cas", "post_cas_a_" + glue::as_glue(loc) + ".rdata"))
  
  exo <- dat[c("r_death_a", "r_death_s", "ppv_pub", "ppv_eng", "sens_acf", "spec_acf", "dur_tx")]
  pp <- rstan::extract(post, pars=names.pars) %>% 
    as_tibble() %>% 
    rename(wt = "wt[1]") %>% 
    mutate(across(c(starts_with("prv_")), function(x) x * wt)) %>% 
    select(-wt) %>% 
    bind_cols(exo)
  
  ps <- lapply(i_sel, function(i) as.list(pp[i, ]))
  ps_a000 <- lapply(ps, transform_by_ratio, rat_dx=0)    %>% 
    bind_rows() %>% mutate(Scenario = "R000", Baseline = "Eq", Loc = loc)
  ps_a025 <- lapply(ps, transform_by_ratio, rat_dx=0.25) %>% 
    bind_rows() %>% mutate(Scenario = "R025", Baseline = "Eq", Loc = loc)
  ps_a050 <- lapply(ps, transform_by_ratio, rat_dx=0.5)  %>% 
    bind_rows() %>% mutate(Scenario = "R050", Baseline = "Eq", Loc = loc)
  ps_a075 <- lapply(ps, transform_by_ratio, rat_dx=0.75) %>% 
    bind_rows() %>% mutate(Scenario = "R075", Baseline = "Eq", Loc = loc)
  ps_a100 <- lapply(ps, transform_by_ratio, rat_dx=1) %>% 
    bind_rows() %>% mutate(Scenario = "R100", Baseline = "Eq", Loc = loc)
  
  
  load(here::here("out", "sub_cas", "post_cas_b_" + glue::as_glue(loc) + ".rdata"))
  
  exo <- dat[c("r_death_a", "r_death_s", "ppv_pub", "ppv_eng", "sens_acf", "spec_acf", "dur_tx")]
  pp <- rstan::extract(post, pars=names.pars) %>% 
    as_tibble() %>% 
    rename(wt = "wt[1]") %>% 
    mutate(across(c(starts_with("prv_")), function(x) x * wt)) %>% 
    select(-wt) %>% 
    bind_cols(exo)
  
  ps <- lapply(i_sel, function(i) as.list(pp[i, ]))
  ps_b000 <- lapply(ps, transform_by_ratio, rat_dx=0)    %>% 
    bind_rows() %>% mutate(Scenario = "R000", Baseline = "Zero0", Loc = loc)
  ps_b025 <- lapply(ps, transform_by_ratio, rat_dx=0.25) %>% 
    bind_rows() %>% mutate(Scenario = "R025", Baseline = "Zero0", Loc = loc)
  ps_b050 <- lapply(ps, transform_by_ratio, rat_dx=0.5)  %>% 
    bind_rows() %>% mutate(Scenario = "R050", Baseline = "Zero0", Loc = loc)
  ps_b075 <- lapply(ps, transform_by_ratio, rat_dx=0.75) %>% 
    bind_rows() %>% mutate(Scenario = "R075", Baseline = "Zero0", Loc = loc)
  ps_b100 <- lapply(ps, transform_by_ratio, rat_dx=1) %>% 
    bind_rows() %>% mutate(Scenario = "R100", Baseline = "Zero0", Loc = loc)
  
  
  load(here::here("out", "sub_cas", "post_cas_c_" + glue::as_glue(loc) + ".rdata"))
  
  exo <- dat[c("r_death_a", "r_death_s", "ppv_pub", "ppv_eng", "sens_acf", "spec_acf", "dur_tx")]
  pp <- rstan::extract(post, pars=names.pars) %>% 
    as_tibble() %>% 
    rename(wt = "wt[1]") %>% 
    mutate(across(c(starts_with("prv_")), function(x) x * wt)) %>% 
    select(-wt) %>% 
    bind_cols(exo)
  
  ps <- lapply(i_sel, function(i) as.list(pp[i, ]))
  ps_c000 <- lapply(ps, transform_by_ratio, rat_dx=0)    %>% 
    bind_rows() %>% mutate(Scenario = "R000", Baseline = "Free", Loc = loc)
  ps_c025 <- lapply(ps, transform_by_ratio, rat_dx=0.25) %>% 
    bind_rows() %>% mutate(Scenario = "R025", Baseline = "Free", Loc = loc)
  ps_c050 <- lapply(ps, transform_by_ratio, rat_dx=0.5)  %>% 
    bind_rows() %>% mutate(Scenario = "R050", Baseline = "Free", Loc = loc)
  ps_c075 <- lapply(ps, transform_by_ratio, rat_dx=0.75) %>% 
    bind_rows() %>% mutate(Scenario = "R075", Baseline = "Free", Loc = loc)
  ps_c100 <- lapply(ps, transform_by_ratio, rat_dx=1) %>% 
    bind_rows() %>% mutate(Scenario = "R100", Baseline = "Free", Loc = loc)
  
  
  
  transformed <- bind_rows(
    ps_a000, ps_a025, ps_a050, ps_a075, ps_a100,
    ps_b000, ps_b025, ps_b050, ps_b075, ps_b100,
    ps_c000, ps_c025, ps_c050, ps_c075, ps_c100
  ) %>% 
    filter(r_recsi > 0)
  
  
  sims <- bind_rows(lapply(1:nrow(transformed), function(i) {
    p <- as.list(transformed[i, ])
    
    cm$set_user(user=p)
    cm$set_user(intv = 0)
    
    ys_0 <- cm$run(2019:2025)
    
    cm$set_user(intv = 0.5)
    
    ys_1 <- cm$run(2019:2025)
    
    res <- list(
      Delay0 = ys_0[nrow(ys_0), 'del_tot'],
      Delay1 = ys_1[nrow(ys_1), 'del_tot'],
      CDR0 = ys_0[nrow(ys_0), 'cdr'],
      CDR1 = ys_1[nrow(ys_1), 'cdr']
    )
    
    c(p, res)
    
  }))
  
  write_csv(sims, here::here("out", "sub_cas", "sims_" + glue::as_glue(loc) + ".csv"))
}



