library(tidyverse)


source(here::here("R", "calc_cascade.R"))


for (scenario in c("shared_pr_asym", "shared_r_onset")) {
  for (cnr_year in 2019:2021) {
    folder <- paste0("cas_", cnr_year)
    file_rates <- glue::as_glue("rates_") + scenario + ".rdata"
    
    load(here::here("out", folder, file_rates))
    
    
    cascades <- rates %>% calc_cascade()
    
    save(cascades, settings, file = here::here("out", folder, glue::as_glue("cascades_") + scenario + ".rdata"))
  }
}
