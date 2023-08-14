library(tidyverse)



cas <- read_csv(here::here("docs", "tabs", "cascade.csv"))


tab_trs <- cas %>% 
  select(Location, starts_with("p_")) %>% 
  group_by(Location) %>% 
  summarise(across(everything(), mean))


tab_trs %>% 
  write_csv(here::here("docs", "tabs", "trs.csv"))



tab_end <- cas %>% 
  mutate(
    EndDie = DropDieA + DropDieS + DropDieC + DropDieT,
    EndLTFU = DropLTFU,
    EndSelfCure = DropSelfCureA + DropSelfCureS + DropSelfCureC,
    EndTxSucc = CasTxS
  ) %>% 
  select(Location, starts_with("End")) %>% 
  group_by(Location) %>% 
  summarise(across(everything(), mean))

tab_end %>% 
  write_csv(here::here("docs", "tabs", "end.csv"))


tab_delay <- cas %>% 
  mutate(
    DelayPat = DelayPat * 12,
    DelaySys = DelaySys * 12,
    DelayTot = DelayTot * 12
  ) %>% 
  select(Location, starts_with("Delay"), starts_with("Cas")) %>% 
  group_by(Location) %>% 
  summarise(across(everything(), mean))

tab_delay %>% 
  write_csv(here::here("docs", "tabs", "delay.csv"))

