library(tidyverse)
library(tidybayes)


source(here::here("R", "calc_shifting.R"))


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

load(here::here("data", "shifting", "shifting_mat.rdata"))


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



fn(c(10, 0.9, 0.5, 1.2), mat0 = mat0, det = det, pdi_pri = 1)



extract_shifting(10, 0.9, 0.5, 1.2, mat0 = mat0, det = det)$dy %>% 
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



xx <- bind_rows(lapply(seq(0.1, 0.95, 0.05), function(pdi_pri) {
  
  opt <- optim(c(0.9, pdi_pri, 0.9, 1), fn, method = "L-BFGS-B", lower = rep(0.1, 4), upper = c(50, 0.99999, 0.99, 20), 
               mat0 = mat0, det = det, pdi_pri = pdi_pri)
  
  x = unname(opt$par)
  list(
    amp_pri=x[1],
    pdi_pri=pdi_pri,
    p_eng = x[3],
    adj_pub = x[4]
  )
}))



plot(xx$amp_pri, xx$pdi_pri)
plot(xx$p_eng, xx$pdi_pri)




opt <- optim(c(0.9, 0.5, 0.9, 1), fn, method = "L-BFGS-B", lower = rep(0.1, 4), upper = c(50, 0.99999, 0.99, 20), 
             mat0 = mat0, det = det, pdi_pri = 0.5)

sim <- extract_shifting(opt$par[1], pdi_pri = 0.5, opt$par[3], opt$par[4], mat0 = mat0, det = det)



sim$dy %>% 
  mutate(
    pub_det = cumsum(pub_det),
    eng_det = cumsum(eng_det),
    pri_det = cumsum(pri_det)
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


save(sim, file = here::here("out", "shifting.rdata"))
