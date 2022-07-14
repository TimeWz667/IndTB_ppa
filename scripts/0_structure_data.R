library(tidyverse)



dat_state <- read_csv(here::here("data", "StateMap.csv")) %>% 
  mutate(
    State = ifelse(State == "Chhatisgarh", "Chhattisgarh", State),
    StateGroup = ifelse(StateGroup == "Chhatisgarh", "Chhattisgarh", StateGroup)
  )


dat_nikshay <- read_csv(here::here("data", "Nikshay.csv")) %>% 
  mutate(
    State = ifelse(State == "Jammu Kashmir", "Jammu & Kashmir", State)
  ) %>% 
  filter(State %in% dat_tbps$State)


dat_pop <- dat_nikshay %>% 
  select(State, Pop) %>% 
  left_join(dat_state %>% select(State, Region))


drug <- read_csv(here::here("data", "drug_sale.csv")) %>% 
  extract(Tx_month_pri, c("TxMo_M", "TxMo_L", "TxMo_U"), "(\\S+) \\((\\S+), (\\S+)\\)", convert = T) %>% 
  mutate(
    across(starts_with("TxMo"), function(x) x * 1e-5 / 12),
    txmo0 = 8.39 * 1e5 * 7 / 1.35e9 / 12 / TxMo_M[State == "India"],
    DrugTime_M = TxMo_M * txmo0, 
    DrugTime_L = TxMo_L * txmo0, 
    DrugTime_U = TxMo_U * txmo0,
    State = case_when(
      State == "Jammu and Kashmir" ~ "Jammu & Kashmir",
      State == "Tamilnadu" ~ "Tamil Nadu",
      T ~ State
    )
  ) %>% 
  full_join(dat_pop) %>% 
  mutate(
    DrugTime_M = ifelse(is.na(DrugTime_M), DrugTime_M[State == "North East"], DrugTime_M),
    DrugTime_L = ifelse(is.na(DrugTime_L), DrugTime_L[State == "North East"], DrugTime_L),
    DrugTime_U = ifelse(is.na(DrugTime_U), DrugTime_U[State == "North East"], DrugTime_U)
  ) %>% 
  filter(State != "North East") %>% 
  select(Region, State, starts_with("DrugTime"))



dat_tbps <- read_csv(here::here("data", "TBPS_State.csv")) %>% 
  mutate(
    State = ifelse(State == "Chhatisgarh", "Chhattisgarh", State),
    StateGroup = ifelse(StateGroup == "Chhatisgarh", "Chhattisgarh", StateGroup)
  )



dat_nikshay2019 <- dat_nikshay %>% 
  select(State, Pop,
         CNR_pub = CNR_public_2019, CNR_eng = CNR_private_2019,
         TxI_pub = TxI_public_2019, TxI_eng = TxI_private_2019,
         TxSucc_pub = TxSucc_public_2019, TxSucc_eng = TxSucc_private_2019,
         TxDead_pub = TxDead_public_2019, TxDead_eng = TxDead_private_2019,
         TxLTFU_pub = TxLTFU_public_2019, TxLTFU_eng = TxLTFU_private_2019) %>% 
  mutate(
    across(c(CNR_pub, CNR_eng), function(x) x * 1e-5),
    across(c(TxI_pub, TxI_eng), function(x) x * 1e-2),
    Tx_pub = TxSucc_pub + TxLTFU_pub + TxDead_pub,
    across(c(TxSucc_pub, TxLTFU_pub, TxDead_pub), function(x) x / Tx_pub),
    Tx_eng = TxSucc_eng + TxLTFU_eng + TxDead_eng,
    across(c(TxSucc_eng, TxLTFU_eng, TxDead_eng), function(x) x / Tx_eng)
  ) %>% 
  select(-Tx_pub, -Tx_eng)


dat_nikshay2020 <- dat_nikshay %>% 
  select(State, Pop,
         CNR_pub = CNR_public_2020, CNR_eng = CNR_private_2020,
         TxI_pub = TxI_public_2019, TxI_eng = TxI_private_2019,
         TxSucc_pub = TxSucc_public_2019, TxSucc_eng = TxSucc_private_2019,
         TxDead_pub = TxDead_public_2019, TxDead_eng = TxDead_private_2019,
         TxLTFU_pub = TxLTFU_public_2019, TxLTFU_eng = TxLTFU_private_2019) %>% 
  mutate(
    across(c(CNR_pub, CNR_eng), function(x) x * 1e-5),
    across(c(TxI_pub, TxI_eng), function(x) x * 1e-2),
    Tx_pub = TxSucc_pub + TxLTFU_pub + TxDead_pub,
    across(c(TxSucc_pub, TxLTFU_pub, TxDead_pub), function(x) x / Tx_pub),
    Tx_eng = TxSucc_eng + TxLTFU_eng + TxDead_eng,
    across(c(TxSucc_eng, TxLTFU_eng, TxDead_eng), function(x) x / Tx_eng)
  ) %>% 
  select(-Tx_pub, -Tx_eng)


dat_nikshay2021 <- dat_nikshay %>% 
  select(State, Pop,
         CNR_pub = CNR_public_2021, CNR_eng = CNR_private_2021,
         TxI_pub = TxI_public_2019, TxI_eng = TxI_private_2019,
         TxSucc_pub = TxSucc_public_2019, TxSucc_eng = TxSucc_private_2019,
         TxDead_pub = TxDead_public_2019, TxDead_eng = TxDead_private_2019,
         TxLTFU_pub = TxLTFU_public_2019, TxLTFU_eng = TxLTFU_private_2019) %>% 
  mutate(
    across(c(CNR_pub, CNR_eng), function(x) x * 1e-5),
    across(c(TxI_pub, TxI_eng), function(x) x * 1e-2),
    Tx_pub = TxSucc_pub + TxLTFU_pub + TxDead_pub,
    across(c(TxSucc_pub, TxLTFU_pub, TxDead_pub), function(x) x / Tx_pub),
    Tx_eng = TxSucc_eng + TxLTFU_eng + TxDead_eng,
    across(c(TxSucc_eng, TxLTFU_eng, TxDead_eng), function(x) x / Tx_eng)
  ) %>% 
  select(-Tx_pub, -Tx_eng)
  

d_tbps <- dat_tbps %>% 
  mutate(
    PrevUt = N_TB0 / N * AmpAll,
    PrevTx = N_OnATT2 / N * AmpAll,
    PrevTxPub = N_OnATT_Pub / N * AmpAll,
    PrevTxPri = N_OnATT_Pri / N * AmpAll,
    Pr_Asym = N_Asym / N_TB0,
    Sym = NotAware0 + NotCS0 + NotDet0,
    Pr_NotAware = (1 - Pr_Asym) * NotAware0 / Sym,
    Pr_NotCS = (1 - Pr_Asym) * NotCS0 / Sym,
    Pr_NotDet = (1 - Pr_Asym) * NotDet0 / Sym
  ) %>% 
  select(Region, State, N, N_TB0, starts_with("Prev"), starts_with("Pr_"))




d_cascade <- d_tbps %>% left_join(dat_nikshay2019) %>% left_join(drug)
save(d_cascade, file = here::here("data", "cascade", "d_cascade_2019.rdata"))
write_csv(d_cascade, here::here("data", "cascade", "d_cascade_2019.csv"))


d_cascade <- d_tbps %>% left_join(dat_nikshay2020) %>% left_join(drug)
save(d_cascade, file = here::here("data", "cascade", "d_cascade_2020.rdata"))
write_csv(d_cascade, here::here("data", "cascade", "d_cascade_2020.csv"))


d_cascade <- d_tbps %>% left_join(dat_nikshay2021) %>% left_join(drug)
save(d_cascade, file = here::here("data", "cascade", "d_cascade_2021.rdata"))
write_csv(d_cascade, here::here("data", "cascade", "d_cascade_2021.csv"))

