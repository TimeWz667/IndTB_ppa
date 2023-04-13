library(tidyverse)



loc <- "India"


pars <- read_csv(file=here::here("out", "sub_cs", "post_cs_s0_" + glue::as_glue(loc) + ".csv"))


re = pars %>% 
  mutate(
    det = det_eng + det_pub + det_pri,
    fp = det_pub * (1 / ppv_pub - 1) + det_eng * (1 / ppv_eng - 1) + det_pri * (1 / ppv_pri - 1),
    p_under = det_pri / det,
    ppv = det / (det + fp)
  )
summary(re)


## Summarise results
tab <- bind_rows(lapply(locations, function(loc) {
  read_csv(file=here::here("out", "sub_cs", "post_cs_s1_" + glue::as_glue(loc) + ".csv")) %>% 
    summarise(across(everything(), list(
      m = mean,
      l = function(x) quantile(x, 0.025),
      u = function(x) quantile(x, 0.975)
    ))) %>% 
    mutate(State = loc) %>%
    relocate(State)
  
}))
write_csv(tab, here::here("docs", "tabs", "post_cs_s1.csv"))


tab <- bind_rows(lapply(locations, function(loc) {
  read_csv(file=here::here("out", "sub_cs", "post_cs_s2_" + glue::as_glue(loc) + ".csv")) %>% 
    summarise(across(everything(), list(
      m = mean,
      l = function(x) quantile(x, 0.025),
      u = function(x) quantile(x, 0.975)
    ))) %>% 
    mutate(State = loc) %>%
    relocate(State)
  
}))
write_csv(tab, here::here("docs", "tabs", "post_cs_s2.csv"))


tab <- bind_rows(lapply(locations, function(loc) {
  read_csv(file=here::here("out", "sub_cs", "post_cs_s3_" + glue::as_glue(loc) + ".csv")) %>% 
    summarise(across(everything(), list(
      m = mean,
      l = function(x) quantile(x, 0.025),
      u = function(x) quantile(x, 0.975)
    ))) %>% 
    mutate(State = loc) %>%
    relocate(State)
  
}))
write_csv(tab, here::here("docs", "tabs", "post_cs_s3.csv"))
