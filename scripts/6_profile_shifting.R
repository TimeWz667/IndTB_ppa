library(tidyverse)
library(odin)


source(here::here("R", "calc_shifting.R"))
model_shifting <- odin::odin("odin/shifting.R", verbose = F) 



for (cnr_year in 2019:2021) {
  folder <- glue::as_glue("cas_") + cnr_year
  load(here::here("out", folder, "shifting.rdata"))

  pars_shifting0 <- lapply(sim_shifting, function(s0) {
    s0 <- s0$Sim
    trm <- s0$tr_mat[, 1:3] + s0$tr_mat[, 4:6]
    
    list(
      entry = trm["UC", ],
      p_diag0 = s0$p_diag["UC", ],
      p_diag1 = s0$p_diag[1:3, ],
      trm_shift = trm[1:3, ]
    )
  })

  
  for (scenario in c("shared_pr_asym", "shared_r_onset")) {
    file_rates <- glue::as_glue("rates_") + scenario + ".rdata"
    file_shift <- glue::as_glue("shifting_") + scenario + ".rdata"
    
    load(here::here("out", folder, file_rates))

    locs <- unique(rates$Location)
    names(locs) <- locs
    
    pars_shifting <- lapply(locs, function(loc) {
      print(loc)
      rates_loc <- rates %>% filter(Location == loc)
      shift_loc <- sim_shifting[[loc]]$Sim
      
      list(
        Rates = bind_rows(lapply(unique(rates_loc$Key), function(k) {
          d <- rates_loc %>% filter(Key == k) %>% as.list()
          d$r_recsi <- find_r_recsi(d, sim_shifting = shift_loc)
          d
        })),
        Shifting = pars_shifting0[[loc]]
      )
    })
    save(pars_shifting, file = here::here("out", folder, file_shift))  
  }
}





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
