library(tidyverse)


m_dy <- odin::odin("odin/m_dy.R")



cnr_year <- 2019
folder <- glue::as_glue("cas_") + cnr_year

scenario <- "shared_pr_asym"
file_shift <- glue::as_glue("shifting_") + scenario + ".rdata"


load(file = here::here("out", folder, file_shift))


pars <- pars_shifting$Andhra_Pradesh

k = 1
d <- pars$Rates %>% filter(Key == k) %>% as.list()


inp <- c(list(
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


cm <- m_dy$new(user = inp)


ys <- cm$run(seq(0, 10, 0.1))

ys[nrow(ys), c("det_pub", "det_eng", "det_pri")]
c(d$DetR_pub, d$DetR_eng, d$DetR_pri)


ys[, "PrevUT"]


