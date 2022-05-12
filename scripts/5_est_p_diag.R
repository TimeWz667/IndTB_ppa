library(tidyverse)


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

n_visits <- 2.35 # from Muniyandi et al. 2020



##
load(here::here("data", "shifting_mat.rdata"))


mat0 <- tapply(tr_mat$Pr, list(tr_mat$From, tr_mat$To), sum)
mat0[is.na(mat0)] <- 0


extract_shifting <- function(x, mat0, det) {
  mat <- mat0
  mat <- cbind(mat, pri_det = 0)
  # mat["pub", "pub_det"] <- mat["pub", "pub_det"] * x[1]
  mat["pub", "pub_det"] <- x[1]
  mat["pub", 1:2] <- mat["pub", 1:2] / sum(mat["pub", 1:2]) * (1 - mat["pub", "pub_det"])
  
  p_diag_pub <- mat['pub', 'pub_det']
  p_diag_pri <- x[2] * p_diag_pub
  
  mat["pri", ] <- mat["pri", ] * (1 - p_diag_pri)
  mat["pri", "pri_det"] <- p_diag_pri
  
  mat["UC", c("pub", "pri")] <- mat["UC", c("pub", "pri")] * c(1, 1 / (1 - p_diag_pri))
  mat["UC", ] <- mat["UC", ] / sum(mat["UC", ])
  
  
  p_diag_pub <- p_diag_eng <- mat['pub', 'pub_det']
  p_diag_pri0 <- mat["pri", "pri_det"]
  
  eng <- det["eng"] / sum(det[2:3]) * p_diag_pri0 / p_diag_eng
  
  p_diag_pri <- (p_diag_pri0 - eng * p_diag_eng)/(1 - eng)
  p_diag_pri <- unname(p_diag_pri)
  

  p_diag <- c(p_diag_pub, p_diag_eng, p_diag_pri)
  
  ent_pri <- det[2:3] / p_diag[2:3]
  ent_pri <- ent_pri / sum(ent_pri)
  
  tr_mat <- matrix(0, 4, 6)
  colnames(tr_mat) <- c("pub", "eng", "pri", "pub_det", "eng_det", "pri_det")
  rownames(tr_mat) <- c("pub", "eng", "pri", "UC")
  tr_mat["UC", "pub"] <- mat["UC", "pub"]
  tr_mat["UC", c("eng", "pri")] <- mat["UC", "pri"] * ent_pri
  tr_mat["UC", ] <- c(tr_mat["UC", 1:3] * (1 - p_diag), tr_mat["UC", 1:3] * p_diag)
  tr_mat["pub", c("pub", "pub_det")] <- mat["pub", c("pub", "pub_det")]
  tr_mat["pub", c("eng", "pri")] <- mat["pub", "pri"] * ent_pri
  tr_mat["pri", c("pub", "pri_det")] <- mat["pri", c("pub", "pri_det")]
  tr_mat["pri", c("eng", "pri")] <- mat["pri", "pri"] * ent_pri
  tr_mat["eng", c("pub", "eng_det")] <- mat["pri", c("pub", "pri_det")]
  tr_mat["eng", c("eng", "pri")] <- mat["pri", "pri"] * ent_pri

  tr_mat["eng", "eng_det"] <- p_diag[2]
  tr_mat["eng", 1:3] <- tr_mat["eng", 1:3] / sum(tr_mat["eng", 1:3]) * (1 - p_diag[2])
  tr_mat["pri", "pri_det"] <- p_diag[3]
  tr_mat["pri", 1:3] <- tr_mat["pri", 1:3] / sum(tr_mat["pri", 1:3]) * (1 - p_diag[3])
  
  m0 <- list(vis = 1)
  m0 <- c(m0, tr_mat["UC", ])
  
  ms <- list(m0)
  
  for (i in 2:20) {
    m <- list(vis = i)
    m <- c(m, colSums(tr_mat[c("pub", "eng", "pri"), ] * c(m0$pub, m0$eng, m0$pri)))

    ms[[length(ms) + 1]] <- m
    m0 <- m
  }
  
  ms <- bind_rows(ms)
  
  return(list(
    tr_mat = tr_mat,
    p_diag = p_diag,
    dy = ms
  ))
  
  return (ms)
}


fn <- function(x, mat0 = mat0, det = det) {
  sim <- extract_shifting(x, mat0 = mat0, det = det)$dy %>% 
    summarise(
      n_vis = sum(vis * pub_det) / sum(pub_det),
      pub = sum(pub_det),
      eng = sum(eng_det),
      pri = sum(pri_det),
      det = pub + eng + pri,
      pub_det = pub / det
    ) %>% 
    select(n_vis, pub_det) %>% 
    unlist()
  
  sum((sim / c(2.35, 0.6307134) - 1) ^ 2)
}



opt <- optim(c(0.9, 0.9), fn, method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(1, 1), mat0 = mat0, det = det)




sim <- extract_shifting(opt$par, mat0 = mat0, det = det)


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
  ggplot() +
  geom_rect(aes(xmin = vis - 0.4, xmax = vis + 0.4, ymin = y0, ymax = y1, fill = name)) +
  theme_bw()


sim$dy %>% 
  summarise(
    n_vis = sum(vis * pub_det) / sum(pub_det),
    pub = sum(pub_det),
    eng = sum(eng_det),
    pri = sum(pri_det),
    det = pub + eng + pri,
    pub_det = pub / det
  ) %>% 
  select(n_vis, pub_det) %>% 
  unlist()


sim$dy %>% 
  mutate(
    pub_det = cumsum(pub_det),
    eng_det = cumsum(eng_det),
    pri_det = cumsum(pri_det)
  ) %>% 
  select(ends_with("_det"))


