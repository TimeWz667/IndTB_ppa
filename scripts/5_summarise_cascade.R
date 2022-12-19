library(tidyverse)



ds <- dir(here::here("out", "sim", "cas"))
ds <- ds[startsWith(ds, "pars_onds")]


lvs <- c(
  "S_Asym", "S_Sym", "S_ExCs", "S_TxPub", "S_TxPri", 
  "S_Cured", "S_SelfCured", "S_LTFU", "S_DeadU", "S_DeadT"
)


file <- "pars_onds_Andhra_Pradesh.csv"


sim <- read.csv(here::here("out", "sim", "cas", file))



sim %>% 
  mutate(Time = Time - min(Time)) %>% 
  group_by(Source, Time) %>% 
  filter(Time <= 2) %>% 
  summarise(across(starts_with("S_"), mean)) %>% 
  ungroup() %>% 
  pivot_longer(starts_with("S_"), names_to = "Stage") %>% 
  mutate(
    Stage = factor(Stage, rev(lvs))
  ) %>% 
  ggplot() +
  geom_bar(aes(x = Time, y = value, fill = Stage), position=position_stack(), stat="identity")



sim %>% 
  group_by(Key) %>% 
  summarise(
    Asym = sum(S_Asym) * 0.1 * 12,
    Sym = sum(S_Sym) * 0.1 * 12,
    ExCs = sum(S_ExCs) * 0.1 * 12,
  ) %>% 
  ungroup() %>% 
  summarise(across(everything(), list(
    M = mean,
    L = function(x) quantile(x, 0.025),
    U = function(x) quantile(x, 0.975)
  )))






