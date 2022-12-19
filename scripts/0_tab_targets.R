library(tidyverse)



load(here::here("data", "dat_cas.rdata"))


targets <- bind_rows(
  dat_noti %>% 
    mutate(
      CNR_Pub = N_Det_Pub / Pop,
      CNR_Eng = N_Det_Eng / Pop,
      CNR_Acf = N_Det_ACF / Pop
    ) %>% 
    select(State, Region, Year, starts_with("CNR_")) %>% 
    pivot_longer(-c(State, Region, Year), names_to = "Index", values_to = "Observed"),
  dat_tbps %>% 
    mutate(
      PrevA = N_Asym / N,
      PrevS = N_NotCS / N,
      PrevC = N_NotDet / N,
      PrevTxPri = N_OnATT_Pri / N,
      OnPriDrugM = DrugTime_M,
      OnPriDrugL = DrugTime_L,
      OnPriDrugU = DrugTime_U,
      Year = 2020
    ) %>% 
    select(State, Region, Year, starts_with("Prev"), starts_with("OnPri")) %>% 
    pivot_longer(-c(State, Region, Year), names_to = "Index", values_to = "Observed")
)


targets %>% write_csv(here::here("data", "targets.csv"))



