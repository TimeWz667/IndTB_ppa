library(tidyverse)


dat_state <- read_csv(here::here("data", "StateMap.csv"))

dat_nikshay <- read_csv(here::here("data", "Nikshay.csv")) %>% 
  mutate(
    State = ifelse(State == "Jammu Kashmir", "Jammu & Kashmir", State)
  )


dat_tbps <- read_csv(here::here("data", "TBPS_State.csv")) %>% 
  mutate(
    State = ifelse(State == "Chhatisgarh", "Chhattisgarh", State),
    StateGroup = ifelse(StateGroup == "Chhatisgarh", "Chhattisgarh", StateGroup)
  )



dat_nikshay2019 <- dat_nikshay %>% 
  filter(State %in% dat_tbps$State) %>% 
  select(State, 
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
  filter(State %in% dat_tbps$State) %>% 
  select(State, 
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
  filter(State %in% dat_tbps$State) %>% 
  select(State, 
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
    Prev = N_TB0 / N * AmpAll,
    PrevTx = N_OnATT2 / N * AmpAll,
    Prev_Std = sqrt(Prev * (1 - Prev) / N),
    PrevTx_Std = sqrt(PrevTx * (1 - PrevTx) / N),
    across(c(Prev, Prev_Std, PrevTx, PrevTx_Std), function(x) x * 1e5),
    Pr_Asym = N_Asym / N_TB0,
    Sym = NotAware0 + NotCS0 + NotDet0,
    Pr_NotAware = (1 - Pr_Asym) * NotAware0 / Sym,
    Pr_NotCS = (1 - Pr_Asym) * NotCS0 / Sym,
    Pr_NotDet0 = (1 - Pr_Asym) * NotDet0 / Sym
  ) %>% 
  select(Region, State, starts_with("Prev"), starts_with("Pr_"))




d_cascade <- d_tbps %>% left_join(dat_nikshay2019)

save(d_cascade, file = here::here("data", "cascade", "d_cascade_2019.rdata"))


d_cascade <- d_tbps %>% left_join(dat_nikshay2020)

save(d_cascade, file = here::here("data", "cascade", "d_cascade_2020.rdata"))


d_cascade <- d_tbps %>% left_join(dat_nikshay2021)

save(d_cascade, file = here::here("data", "cascade", "d_cascade_2021.rdata"))

