library(tidyverse)
library(tidybayes)


source(here::here("R", "calc_shifting.R"))


assumed_pdi_pri <- 0.5

## Prepare data target
load(here::here("out", "cas_2019", "rates_shared_pr_asym.rdata"))


det <- rates %>% 
  filter(Location == "Tamil_Nadu") %>% 
  select(starts_with("DetR_")) %>% 
  distinct() %>% 
  pivot_longer(everything()) %>% 
  separate(name, c("Index", "sector"), "_") %>% 
  mutate(
    det = value / sum(value)
  )
det <- set_names(det$det, det$sector)
det

n_visits <- 2.75 # from Muniyandi et al. 2020



##

load(here::here("data", "shifting_mat.rdata"))


mat0 <- tapply(tr_mat$Pr, list(tr_mat$From, tr_mat$To), sum)
mat0[is.na(mat0)] <- 0


fn <- function(x, mat0 = mat0, det = det, pdi_pri=NA) {
  amp_pri <- x[1]
  pdi_pri <- ifelse(is.na(pdi_pri), x[2], pdi_pri)
  p_eng <- x[3]
  adj_pub <- x[4]
  
  sim <- extract_shifting(amp_pri, pdi_pri, p_eng, adj_pub, mat0 = mat0, det = det)$dy %>% 
    summarise(
      n_vis = sum(vis * pub_det) / sum(pub_det),
      pub = sum(pub_det),
      eng = sum(eng_det),
      pri = sum(pri_det),
      det = pub + eng + pri,
      pub_det = pub / det,
      eng_det = eng / det
    ) %>% 
    select(n_vis, pub_det, eng_det) %>% 
    unlist()
  
  sum((sim / c(n_visits, det[1:2]) - 1) ^ 2)
}


opt <- optim(c(0.9, assumed_pdi_pri, 0.9, 1), fn, method = "L-BFGS-B", 
             lower = rep(0.1, 4), upper = c(50, 0.99999, 0.99, 20), 
             mat0 = mat0, det = det, pdi_pri = assumed_pdi_pri)

sim <- extract_shifting(opt$par[1], pdi_pri = assumed_pdi_pri, 
                        opt$par[3], opt$par[4], mat0 = mat0, det = det)

sol <- opt$par
sol[2] <- assumed_pdi_pri
names(sol) <- c("amp_pri", "pdi_pri", "p_eng", "adj_pub")


sim$dy %>% 
  mutate(
    across(c(pub_det, eng_det, pri_det), cumsum)
  ) %>% 
  pivot_longer(-vis) %>% 
  group_by(vis) %>% 
  mutate(
    y1 = cumsum(value),
    y0 = c(0, y1[-length(y1)])
  ) %>% 
  filter(vis < 7) %>% 
  ggplot() +
  geom_rect(aes(xmin = vis - 0.5, xmax = vis + 0.5, 
                ymin = y0, ymax = y1, fill = name), colour = "grey") +
  theme_bw()


sim$dy %>% 
  summarise(
    n_vis = sum(vis * pub_det) / sum(pub_det),
    pub = sum(pub_det),
    eng = sum(eng_det),
    pri = sum(pri_det),
    det = pub + eng + pri,
    pub_det = pub / det,
    eng_det = eng / det,
    pri_det = pri / det
  ) %>% 
  select(n_vis, pub_det, eng_det, pri_det) %>% 
  unlist()


sim$p_diag
sim$tr_mat

sim_shifting <- sim
save(sim_shifting, file = here::here("out", "shifting.rdata"))


### Apply to all locations
fn <- function(x, mat0 = mat0, det = det, sol = sol) {
  sim <- extract_shifting(
    amp_pri = x[1], 
    pdi_pri = sol["pdi_pri"], 
    p_eng = x[2], 
    adj_pub = sol["adj_pub"], 
    mat0 = mat0, det = det
  )$dy %>% 
    summarise(
      pub = sum(pub_det),
      eng = sum(eng_det),
      pri = sum(pri_det),
      det = pub + eng + pri,
      pub_det = pub / det,
      eng_det = eng / det
    ) %>% 
    select(pub_det, eng_det) %>% 
    unlist()
  
  sum((sim / c(det[1:2]) - 1) ^ 2)
}

locs <- rates$Location %>% unique()
names(locs) <- locs


for (cnr_year in 2019:2021) {
  folder <- paste0("cas_", cnr_year)
  file_rates <- glue::as_glue("rates_") + scenario + ".rdata"
  load(here::here("out", folder, file_rates))
  
  sim_shifting <- lapply(locs, function(loc) {
    det <- rates %>% 
      filter(Location == loc) %>% 
      select(starts_with("DetR_")) %>% 
      distinct() %>% 
      pivot_longer(everything()) %>% 
      separate(name, c("Index", "sector"), "_") %>% 
      mutate(
        det = value / sum(value)
      )
    det <- set_names(det$det, det$sector)
    
    opt <- optim(c(0.9, 0.9), fn, method = "L-BFGS-B", 
                 lower = rep(0.1, 2), upper = c(50, 0.99), 
                 mat0 = mat0, det = det, sol = sol)
    
    
    sim <- extract_shifting(
      amp_pri = opt$par[1], pdi_pri = sol["pdi_pri"], 
      p_eng = opt$par[2], adj_pub = sol["adj_pub"], 
      mat0 = mat0, det = det
    )
    
    stat_sim <- sim$dy %>% 
      summarise(
        pub = sum(pub_det), eng = sum(eng_det), pri = sum(pri_det),
        dets = pub + eng + pri,
        pub_det = pub / dets, eng_det = eng / dets, pri_det = pri / dets
      ) %>% 
      mutate(
        Location = loc,
        det_pub_data = det["pub"],
        det_eng_data = det["eng"],
        det_pri_data = det["pri"],
        det_pub_sim = pub_det,
        det_eng_sim = eng_det,
        det_pri_sim = pri_det,
      ) %>% 
      select(Location, starts_with("det_"))
    
    
    list(Sim = sim, Stat = stat_sim)
  })
  
  save(sim_shifting, file = here::here("out", folder, "shifting.rdata"))
}
