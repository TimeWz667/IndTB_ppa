library(tidyverse)
library(odin)


load(here::here("out", "shifting.rdata"))


model_shifting <- odin::odin("odin/shifting.R", verbose = F) 

# for (scenario in c("shared_pr_asym", "shared_r_onset")) {
#   for (cnr_year in 2019:2021) {

scenario = "shared_pr_asym"
cnr_year = 2019

folder <- glue::as_glue("cas_") + cnr_year

file_rates <- glue::as_glue("rates_") + scenario + ".rdata"

load(here::here("out", folder, file_rates))



d <- rates %>% 
  filter(Location == "Tamil_Nadu") %>% 
  filter(Key == 1) %>% 
  as.list()

sim_shifting <- sim





r_recsi <- find_r_recsi(d, sim_shifting = sim_shifting)





## Check results
pdi <- sim_shifting$p_diag
shift <- sim_shifting$tr_mat[, 1:3] + sim_shifting$tr_mat[, 4:6]

cm <- model_shifting$new(user = list(
  fn0 = fn0 <- d$r_cs * shift["UC", ] * (1 - pdi["UC", ]) * d$Prev_C,
  trm = shift[1:3, ],
  p_diag1 = pdi[1:3, ],
  mu = d$r_mu_sym,
  r_recsi = r_recsi
))


ys <- cm$run(seq(0, 100, 0.1))

ys[nrow(ys), c("det_pub", "det_eng", "det_pri")]; 
det1
ys[nrow(ys), "prev_e"]; 
d$Prev_E    
