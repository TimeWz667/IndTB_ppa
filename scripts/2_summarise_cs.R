
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
