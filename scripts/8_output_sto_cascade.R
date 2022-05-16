library(tidyverse)



stats <- list(
  M = median,
  L = function(x) quantile(x, 0.25),
  U = function(x) quantile(x, 0.75)
)


for(cnr_year in 2019:2021) {
  folder <- glue::as_glue("cas_") + cnr_year
  
  for(scenario in c("shared_pr_asym", "shared_r_onset")) {
    load(file = here::here("out", folder, "sims_pp_" + glue::as_glue(scenario) + ".rdata"))

    sel <- sims_pp %>% 
      filter(!is.na(SectorStart) & !is.na(SectorEnd)) %>% 
      mutate(DelPat = DelayPatient * 12, DelSys = DelaySystem * 12)
    
    
    delays <- bind_rows(
      sel %>% 
        group_by(Location, SectorStart) %>% 
        summarise(across(c(DelPat, DelSys), stats)) %>% 
        ungroup() %>% 
        mutate(SectorEnd = "All"),
      sel %>% 
        group_by(Location, SectorEnd) %>% 
        summarise(across(c(DelPat, DelSys), stats)) %>% 
        ungroup() %>% 
        mutate(SectorStart = "All"),
      sel %>% 
        group_by(Location) %>% 
        summarise(across(c(DelPat, DelSys), stats)) %>% 
        ungroup() %>% 
        mutate(SectorEnd = "All", SectorStart = "All"),
      sel %>% 
        group_by(Location, SectorStart, SectorEnd) %>% 
        summarise(across(c(DelPat, DelSys), stats)) %>% 
        ungroup()
    ) %>% 
      mutate(
        SectorStart = factor(SectorStart, c("pub", "eng", "pri", "All")),
        SectorEnd = factor(SectorEnd, c("pub", "eng", "pri", "All")),
        DelPat = sprintf("%.1f (%.1f - %.1f)", DelPat_M, DelPat_L, DelPat_U),
        DelSys = sprintf("%.1f (%.1f - %.1f)", DelSys_M, DelSys_L, DelSys_U)
      ) %>% 
      relocate(Location, SectorStart, SectorEnd)
    
    
    delays %>% 
      filter(Location == "India") %>% 
      (function(ss) {
        tapply(ss$DelSys, list(ss$SectorStart, ss$SectorEnd), function(x) x[1])
      })() %>% 
      write.csv(here::here("results", folder + "_" + scenario, "DelaySto.csv"))
    
    
    
  }
}





