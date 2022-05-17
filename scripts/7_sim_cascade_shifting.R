library(tidyverse)


source(here::here("R", "calc_cascade.R"))


m_cascade <- odin::odin("odin/m_cascade_breakdown.R")


 
for(cnr_year in 2019:2021) {
  folder <- glue::as_glue("cas_") + cnr_year
  print(folder)
  for(scenario in c("shared_pr_asym", "shared_r_onset")) {
    file_shift <- glue::as_glue("shifting_") + scenario + ".rdata"
    load(file = here::here("out", folder, file_shift))
    locs <- names(pars_shifting)

    print(scenario)
    
    cascades_shifting <- bind_rows(lapply(locs, function(loc) {
      print(loc)
      pars <- pars_shifting[[loc]]

      calc_cascade_shifting(pars, m_cascade) %>% mutate(Location = loc)
    }))

    save(cascades_shifting, file = here::here("out", folder, glue::as_glue("cascades_shifting_") + scenario + ".rdata"))
  }
}
