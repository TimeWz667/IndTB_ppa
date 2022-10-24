library(tidyverse)



load(here::here("data", "dat_cas.rdata"))


locations <- unique(dat_noti$State)


locations



tab <- bind_rows(lapply(locations, function(loc) {
  read_csv(file=here::here("out", "sub_cs", "post_cs_nods_" + glue::as_glue(loc) + ".csv")) %>% 
    summarise(across(everything(), list(
      m = mean,
      l = function(x) quantile(x, 0.025),
      u = function(x) quantile(x, 0.975)
    ))) %>% 
    mutate(State = loc) %>%
    relocate(State)
  
}))


write_csv(tab, here::here("docs", "tabs", "post_cs_nods.csv"))



tab <- bind_rows(lapply(locations, function(loc) {
  read_csv(file=here::here("out", "sub_cs", "post_cs_onds_" + glue::as_glue(loc) + ".csv")) %>% 
    summarise(across(everything(), list(
      m = mean,
      l = function(x) quantile(x, 0.025),
      u = function(x) quantile(x, 0.975)
    ))) %>% 
    mutate(State = loc) %>%
    relocate(State)
  
}))


write_csv(tab, here::here("docs", "tabs", "post_cs_onds.csv"))




