library(tidyverse)


d_acf <- local({
  load(here::here("data", "d_itr.rdata"))
  
  d_itr_acf %>% filter(Year == 2020)
})


d_prev <- read_csv(here::here("data", "TBPS", "TBPS_ASC_State.csv"))


d_prev %>% 
  mutate(
    pTBLike = N_TBLike / N,
    pAsym = N_Asym / N,
    pSym = N_NotCS / N,
    pExCs = N_NotDet / N
  ) %>% 
  left_join(d_acf %>% select(State, N_ACF_Mapped)) %>% 
  mutate(
    pVul = round(pmin(N_ACF_Mapped / Pop, 1), 3)
  ) %>% 
  select(State, StateGroup, Region, Pop, pTBLike, pAsym, pSym, pExCs, pVul) %>% 
  write_csv(here::here("data/baseline2020.csv"))



