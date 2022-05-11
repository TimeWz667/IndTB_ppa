library(tidyverse)



data_cascade <- read_csv(here::here("data", "TBPS.csv"))


# With 2019 CNR
d_cascade <- data_cascade %>% 
  filter(!is.na(CNR_public_2021)) %>%
  mutate(
    Pr_NA_NC = NotRecognised + Ignored, 
    Prev_Std = (Prev_U - Prev_L) / 2 / 1.96
  ) %>% 
  select(Location, Pop, Prev = Prev_M, Prev_Std,
         Pr_CS_Sym = Consulted, Pr_NA_NC,
         CNR_pub = CNR_public_2019, CNR_eng = CNR_private_2019,
         TxI_pub = TxI_public_2019, TxI_eng = TxI_private_2019,
         TxSucc_pub = TxSucc_public_2019, TxSucc_eng = TxSucc_private_2019,
         TxDead_pub = TxDead_public_2019, TxDead_eng = TxDead_private_2019,
         TxLTFU_pub = TxLTFU_public_2019, TxLTFU_eng = TxLTFU_private_2019, TxPriPub) %>% 
  mutate(
    across(c(Prev, Prev_Std, CNR_pub, CNR_eng), function(x) x * 1e-5),
    across(c(TxI_pub, TxI_eng), function(x) x * 1e-2),
    Pr_Asym = 0.39,
    Tx_pub = TxSucc_pub + TxLTFU_pub + TxDead_pub,
    across(c(TxSucc_pub, TxLTFU_pub, TxDead_pub), function(x) x / Tx_pub),
    Tx_eng = TxSucc_eng + TxLTFU_eng + TxDead_eng,
    across(c(TxSucc_eng, TxLTFU_eng, TxDead_eng), function(x) x / Tx_eng)
  ) %>% 
  select(-Tx_pub, -Tx_eng)


save(d_cascade, file = here::here("data", "d_cascade_2019.rdata"))



# With 2020 CNR
d_cascade <- data_cascade %>% 
  filter(!is.na(CNR_public_2021)) %>%
  mutate(
    Pr_NA_NC = NotRecognised + Ignored, 
    Prev_Std = (Prev_U - Prev_L) / 2 / 1.96
  ) %>% 
  select(Location, Pop, Prev = Prev_M, Prev_Std,
         Pr_CS_Sym = Consulted, Pr_NA_NC,
         CNR_pub = CNR_public_2020, CNR_eng = CNR_private_2020,
         TxI_pub = TxI_public_2019, TxI_eng = TxI_private_2019,
         TxSucc_pub = TxSucc_public_2019, TxSucc_eng = TxSucc_private_2019,
         TxDead_pub = TxDead_public_2019, TxDead_eng = TxDead_private_2019,
         TxLTFU_pub = TxLTFU_public_2019, TxLTFU_eng = TxLTFU_private_2019, TxPriPub) %>% 
  mutate(
    across(c(Prev, Prev_Std, CNR_pub, CNR_eng), function(x) x * 1e-5),
    across(c(TxI_pub, TxI_eng), function(x) x * 1e-2),
    Pr_Asym = 0.39,
    Tx_pub = TxSucc_pub + TxLTFU_pub + TxDead_pub,
    across(c(TxSucc_pub, TxLTFU_pub, TxDead_pub), function(x) x / Tx_pub),
    Tx_eng = TxSucc_eng + TxLTFU_eng + TxDead_eng,
    across(c(TxSucc_eng, TxLTFU_eng, TxDead_eng), function(x) x / Tx_eng)
  ) %>% 
  select(-Tx_pub, -Tx_eng)


save(d_cascade, file = here::here("data", "d_cascade_2020.rdata"))



# With 2021 CNR
d_cascade <- data_cascade %>% 
  filter(!is.na(CNR_public_2021)) %>%
  mutate(
    Pr_NA_NC = NotRecognised + Ignored, 
    Prev_Std = (Prev_U - Prev_L) / 2 / 1.96
  ) %>% 
  select(Location, Pop, Prev = Prev_M, Prev_Std,
         Pr_CS_Sym = Consulted, Pr_NA_NC,
         CNR_pub = CNR_public_2021, CNR_eng = CNR_private_2021,
         TxI_pub = TxI_public_2019, TxI_eng = TxI_private_2019,
         TxSucc_pub = TxSucc_public_2019, TxSucc_eng = TxSucc_private_2019,
         TxDead_pub = TxDead_public_2019, TxDead_eng = TxDead_private_2019,
         TxLTFU_pub = TxLTFU_public_2019, TxLTFU_eng = TxLTFU_private_2019, TxPriPub) %>% 
  mutate(
    across(c(Prev, Prev_Std, CNR_pub, CNR_eng), function(x) x * 1e-5),
    across(c(TxI_pub, TxI_eng), function(x) x * 1e-2),
    Pr_Asym = 0.39,
    Tx_pub = TxSucc_pub + TxLTFU_pub + TxDead_pub,
    across(c(TxSucc_pub, TxLTFU_pub, TxDead_pub), function(x) x / Tx_pub),
    Tx_eng = TxSucc_eng + TxLTFU_eng + TxDead_eng,
    across(c(TxSucc_eng, TxLTFU_eng, TxDead_eng), function(x) x / Tx_eng)
  ) %>% 
  select(-Tx_pub, -Tx_eng)


save(d_cascade, file = here::here("data", "d_cascade_2021.rdata"))


